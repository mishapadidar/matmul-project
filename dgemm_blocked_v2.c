#include <stdlib.h>
const char* dgemm_desc = "Simple blocked dgemm.";

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 128)
#endif
#ifndef SLICE
#define SLICE ((int) 4)
#endif

/*
  A is M-by-K
  B is K-by-N
  C is M-by-N

  lda is the leading dimension of the matrix (the M of square_dgemm).
*/
void basic_dgemm(const int lda, const int M, const int N, const int K,
                 const double *A, const double *B, double *C)
{
    int i, j, k, l;
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {

            // storage for vector ops
            double *s = (double *) malloc(sizeof(double) * SLICE);
            for (int ndx = 0; ndx < SLICE; ndx++) {
                s[ndx] = 0.0;
            }

            for (k = 0; k < K-(SLICE-1); k+=SLICE) {
                for (l=0; l < SLICE; ++l) {
                    s[l] += A[(k+l)*lda+i] * B[j*lda+k+l];
                }
            }

            double cij = C[j*lda+i];
            // for (k = 0; k < K; ++k) {
            //     cij += A[k*lda+i] * B[j*lda+k];
            // }
            for (int ndx = 0; ndx < SLICE; ndx++)
                cij += s[ndx];
            for (; k < K; ++k)
                cij += A[k*lda+i] * B[j*lda+k];
            C[j*lda+i] = cij;
        }
    }
}

void do_block(const int lda,
              const double *A, const double *B, double *C,
              const int i, const int j, const int k)
{
    const int M = (i+BLOCK_SIZE > lda? lda-i : BLOCK_SIZE);
    const int N = (j+BLOCK_SIZE > lda? lda-j : BLOCK_SIZE);
    const int K = (k+BLOCK_SIZE > lda? lda-k : BLOCK_SIZE);
    basic_dgemm(lda, M, N, K,
                A + i + k*lda, B + k + j*lda, C + i + j*lda);
}

void square_dgemm(const int M, const double *A, const double *B, double *C)
{
    const int n_blocks = M / BLOCK_SIZE + (M%BLOCK_SIZE? 1 : 0);
    int bi, bj, bk;
    for (bi = 0; bi < n_blocks; ++bi) {
        const int i = bi * BLOCK_SIZE;
        for (bj = 0; bj < n_blocks; ++bj) {
            const int j = bj * BLOCK_SIZE;
            for (bk = 0; bk < n_blocks; ++bk) {
                const int k = bk * BLOCK_SIZE;
                do_block(M, A, B, C, i, j, k);
            }
        }
    }
}

