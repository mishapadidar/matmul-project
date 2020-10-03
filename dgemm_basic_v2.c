/*
breaking up loops for vectorization
*/
const char* dgemm_desc = "Basic, three-loop dgemm.";

void square_dgemm(const int M, 
                  const double *A, const double *B, double *C)
{
    int i, j, k, l;
    for (i = 0; i < M; ++i) {
        for (j = 0; j < M; ++j) {
            // double cij = C[j*M+i];
            double s[4] = {0, 0, 0, 0};
            for (k = 0; k < M-3; k+=4) {
                for (l=0; l < 4; ++l) {
                    s[l] += A[(k+l)*M+i] * B[j*M+k+l];
                }
            }
            // Combine sub-sums and handle trailing elements
            double cij = (s[0]+s[1]) + (s[2]+s[3]);
            for (; k < M; ++k)
                cij += A[k*M+i] * B[j*M+k];
            C[j*M+i] = cij;
        }
    }
}
