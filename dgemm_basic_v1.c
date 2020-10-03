/*
Version 1: We switched the loop ordering from column major to row major. We assume that A and B are stored in row-major
format.

Notes: The original version accesses the first matrix in a column-major format and the second matrix in a row-major
format. In order to take advantage of C's proclivity toward's row-major access, we perform the basic
matrix-multiplication accessing the row's of every matrix.
*/
const char* dgemm_desc = "Basic, three-loop dgemm.";

void square_dgemm(const int M, 
                  const double *A, const double *B, double *C)
{
    int i, j, k;
    for (i = 0; i < M; ++i) {
        for (j = 0; j < M; ++j) {
            for (k = 0; k < M; ++k)
                C[i*M+k] += A[j*M+k] * B[i*M+j];
        }
    }
}
