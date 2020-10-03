/*
breaking up loops for vectorization
*/
#include <stdlib.h>

const char* dgemm_desc = "Basic, three-loop dgemm.";

#ifndef SLICE
#define SLICE ((int) 8)
#endif

void square_dgemm(const int M, 
                  const double *A, const double *B, double *C)
{
    int i, j, k, l;
    for (i = 0; i < M; ++i) {
        for (j = 0; j < M; ++j) {
            // double cij = C[j*M+i];
            double *s = (double *) malloc(sizeof(double) * SLICE);
            for (int ndx = 0; ndx < SLICE; ndx++) {
                s[ndx] = 0.0;
            }

            for (k = 0; k < M-(SLICE-1); k+=SLICE) {
                for (l=0; l < SLICE; ++l) {
                    s[l] += A[(k+l)*M+i] * B[j*M+k+l];
                }
            }
            // Combine sub-sums and handle trailing elements
            double cij = 0.0;
            for (int ndx = 0; ndx < SLICE; ndx++)
                cij += s[ndx];
            for (; k < M; ++k)
                cij += A[k*M+i] * B[j*M+k];
            C[j*M+i] = cij;
        }
    }
}
