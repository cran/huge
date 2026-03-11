// blas_config.h — Thin BLAS abstraction for hot-path linear algebra.
//
// When built as an R package, R provides BLAS via $(BLAS_LIBS).
// When standalone, link against system BLAS (OpenBLAS, Accelerate, etc).
// Fortran-style BLAS symbols use trailing underscore on most platforms.
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// y := alpha * A * x + beta * y  (column-major, no-transpose)
// DGEMV: op(A) is m-by-n, x has length n, y has length m
void dgemv_(const char* trans, const int* m, const int* n,
            const double* alpha, const double* A, const int* lda,
            const double* x, const int* incx,
            const double* beta, double* y, const int* incy);

// dot product: sum(x[i]*y[i])
double ddot_(const int* n, const double* x, const int* incx,
             const double* y, const int* incy);

// y := alpha * x + y
void daxpy_(const int* n, const double* alpha,
            const double* x, const int* incx,
            double* y, const int* incy);

// x := alpha * x
void dscal_(const int* n, const double* alpha,
            double* x, const int* incx);

#ifdef __cplusplus
}
#endif
