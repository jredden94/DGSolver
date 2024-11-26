#pragma once

namespace blas {
    void gemv(const double *mat, const double *vec_x, double *vec_y, 
            const unsigned long m, const unsigned long n, const double alpha = 1, const double beta = 0, const bool scale = false);

    void gemm(const double* mat_a, const double *mat_b, double *mat_c, 
            const unsigned long m);

    void InverseMatrix(const double *mat, double *mat_inv, const unsigned short m);

    void Transpose(const double *mat, double *mat_T, const unsigned short m);

    void Transpose(double *mat, const unsigned short n);
}
