// huge_core.h — Standalone C++ core for huge graph estimation algorithms.
// No dependency on Rcpp, pybind11, or Eigen. Pure standard C++.
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace huge {

// ---- Column-major matrix ------------------------------------------------

struct Matrix {
    std::vector<double> v;
    int rows = 0, cols = 0;

    Matrix() = default;
    Matrix(int r, int c) : v(static_cast<size_t>(r) * c, 0.0), rows(r), cols(c) {}

    void resize(int r, int c) { rows = r; cols = c; v.assign(static_cast<size_t>(r) * c, 0.0); }
    void set_zero() { std::fill(v.begin(), v.end(), 0.0); }

    double& operator()(int r, int c)             { return v[static_cast<size_t>(c) * rows + r]; }
    const double& operator()(int r, int c) const { return v[static_cast<size_t>(c) * rows + r]; }
    double* col_ptr(int c)             { return v.data() + static_cast<size_t>(c) * rows; }
    const double* col_ptr(int c) const { return v.data() + static_cast<size_t>(c) * rows; }
};

// ---- MB / TIGER per-column results --------------------------------------

struct ColResult {
    std::vector<double> vals;
    std::vector<int> indices; // encoded as lambda_idx * d + var_idx
};

// ---- Glasso -------------------------------------------------------------

struct GlassoResult {
    std::vector<double> loglik;   // nlambda
    std::vector<double> sparsity; // nlambda
    std::vector<int>    df;       // nlambda
    std::vector<Matrix> path;     // nlambda x (d,d)
    std::vector<Matrix> icov;     // nlambda x (d,d)
    std::vector<Matrix> cov;      // nlambda x (d,d); empty if !cov_output
};

GlassoResult glasso(const double* S_colmajor, int d,
                    const double* lambda, int nlambda,
                    bool scr, bool cov_output);

// ---- MB graph -----------------------------------------------------------

struct MBResult {
    std::vector<ColResult> columns; // size d
};

MBResult mb(const double* S_colmajor, int d,
            const double* lambda, int nlambda);

MBResult mb_scr(const double* S_colmajor, int d,
                const double* lambda, int nlambda,
                const int* idx_scr, int nscr);

// ---- TIGER (sqrt-lasso) -------------------------------------------------

struct TigerResult {
    std::vector<ColResult> columns; // size d
    std::vector<Matrix>    icov;    // nlambda x (d,d)
};

TigerResult tiger(const double* data_colmajor, int n, int d,
                  const double* lambda, int nlambda);

// ---- RIC ----------------------------------------------------------------

double ric(const double* X_colmajor, int n, int d,
           const int* r, int t);

// ---- Scale-free graph generator -----------------------------------------
// G_out: pre-allocated d*d array (column-major), written as adjacency matrix.
// rands: array of (d - d0) uniform(0,1) random values, supplied by caller.

void sfgen(int d0, int d, int* G_out, const double* rands);

} // namespace huge
