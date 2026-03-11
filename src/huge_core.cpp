// huge_core.cpp — Standalone core implementations.
// No Rcpp, no pybind11, no Eigen. Uses BLAS for hot-path linear algebra.
#include "huge/huge_core.h"
#include "huge/blas_config.h"
#include <cstring>

// BLAS constants reused throughout
static const char   BLAS_N   = 'N';
static const int    BLAS_1   = 1;
static const double BLAS_ONE = 1.0, BLAS_ZERO = 0.0;

namespace huge {

static inline double threshold_l1(double x, double thr) {
    if (x > thr) return x - thr;
    if (x < -thr) return x + thr;
    return 0.0;
}

// Column-major accessor for raw const pointer
static inline double cm(const double* data, int nrows, int r, int c) {
    return data[static_cast<size_t>(c) * nrows + r];
}

// =========================================================================
// Glasso: inner coordinate-descent solver
// =========================================================================

static void glasso_sub(Matrix& S, Matrix& W, Matrix& T, int d,
                       double ilambda, int& df, bool scr)
{
    const double thol_act = 1e-4;
    const double thol_ext = 1e-4;
    const int MAX_ITER_EXT = 100;
    const int MAX_ITER_INT = 10000;
    const int MAX_ITER_ACT = 10000;
    const double neg_ilambda = -ilambda;

    // idx_a, idx_i stored as d columns of max-d rows each
    std::vector<int> idx_a(static_cast<size_t>(d) * d, 0);
    std::vector<int> idx_i(static_cast<size_t>(d) * d, 0);
    std::vector<int> size_a(static_cast<size_t>(d), 0);
    std::vector<double> w1(static_cast<size_t>(d), 0.0);
    std::vector<double> ww(static_cast<size_t>(d), 0.0);

    auto idx_a_col = [&](int col) -> int* { return idx_a.data() + static_cast<size_t>(col) * d; };
    auto idx_i_col = [&](int col) -> int* { return idx_i.data() + static_cast<size_t>(col) * d; };

    // Recover initial solution for each individual lasso
    for (int i = 0; i < d; i++) {
        int* ia = idx_a_col(i);
        int* ii = idx_i_col(i);
        double* T_col = T.col_ptr(i);
        const double* S_col = S.col_ptr(i);

        W(i, i) = S_col[i] + ilambda;
        size_a[i] = 0;
        double tmp1 = T_col[i];
        T_col[i] = 0;

        for (int j = 0; j < d; j++) {
            if (scr && std::fabs(S_col[j]) <= ilambda) {
                ii[j] = -1;
                T_col[j] = 0;
                continue;
            }
            if (T_col[j] != 0) {
                ia[size_a[i]++] = j;
                ii[j] = -1;
                T_col[j] = -T_col[j] / tmp1;
            } else {
                ii[j] = 1;
            }
        }
        ii[i] = -1;
    }

    double gap_ext = 1;
    int iter_ext = 0;
    double tmp4 = 0, tmp5 = 0, tmp6 = 0;
    while (gap_ext > thol_ext && iter_ext < MAX_ITER_EXT) {
        tmp6 = 0;
        tmp5 = 0;
        for (int i = 0; i < d; i++) {
            int active_size = size_a[i];
            int* ia = idx_a_col(i);
            int* ii = idx_i_col(i);
            double* T_col = T.col_ptr(i);
            const double* S_col = S.col_ptr(i);

            int gap_int = 1;
            int iter_int = 0;
            for (int j = 0; j < d; j++) ww[j] = T_col[j];

            while (gap_int != 0 && iter_int < MAX_ITER_INT) {
                int size_a_prev = active_size;
                for (int j = 0; j < d; j++) {
                    if (ii[j] != -1) {
                        double r = S_col[j];
                        for (int k = 0; k < active_size; k++)
                            r -= W(ia[k], j) * T_col[ia[k]];

                        double w_new = 0.0;
                        if (r > ilambda) {
                            w_new = (r - ilambda) / W(j, j);
                            ia[active_size++] = j;
                            ii[j] = -1;
                        } else if (r < neg_ilambda) {
                            w_new = (r + ilambda) / W(j, j);
                            ia[active_size++] = j;
                            ii[j] = -1;
                        }
                        w1[j] = w_new;
                        T_col[j] = w_new;
                    }
                }
                gap_int = active_size - size_a_prev;

                double gap_act = 1;
                int iter_act = 0;
                while (gap_act > thol_act && iter_act < MAX_ITER_ACT) {
                    double tmp3 = 0;
                    tmp4 = 0;
                    for (int j = 0; j < active_size; j++) {
                        int w_idx = ia[j];
                        double r = S_col[w_idx] + T_col[w_idx] * W(w_idx, w_idx);
                        for (int k = 0; k < active_size; k++)
                            r -= W(ia[k], w_idx) * T_col[ia[k]];

                        double w_new = 0.0;
                        if (r > ilambda)
                            w_new = (r - ilambda) / W(w_idx, w_idx);
                        else if (r < neg_ilambda)
                            w_new = (r + ilambda) / W(w_idx, w_idx);

                        tmp4 += std::fabs(w_new);
                        tmp3 += std::fabs(w_new - T_col[w_idx]);
                        w1[w_idx] = w_new;
                        T_col[w_idx] = w_new;
                    }
                    gap_act = tmp4 > 0 ? tmp3 / tmp4 : 0;
                    iter_act++;
                }

                // Move false active variables back
                int junk_a = 0;
                for (int j = 0; j < active_size; j++) {
                    int w_idx = ia[j];
                    if (w1[w_idx] == 0) {
                        junk_a++;
                        ii[w_idx] = 1;
                    } else {
                        ia[j - junk_a] = w_idx;
                    }
                }
                active_size -= junk_a;
                iter_int++;
            }
            size_a[i] = active_size;

            // Update W from current T column: W[:,i] = W * T[:,i]
            std::vector<double> temp(d, 0.0);
            dgemv_(&BLAS_N, &d, &d, &BLAS_ONE, W.v.data(), &d,
                   T_col, &BLAS_1, &BLAS_ZERO, temp.data(), &BLAS_1);
            for (int j = 0; j < d; j++) {
                if (j != i) { W(j, i) = temp[j]; W(i, j) = temp[j]; }
            }

            for (int j = 0; j < d; j++)
                tmp5 += std::fabs(ww[j] - T_col[j]);
            tmp6 += tmp4;
        }
        gap_ext = tmp6 > 0 ? tmp5 / tmp6 : 0;
        iter_ext++;
    }

    // Compute final T (precision matrix)
    for (int i = 0; i < d; i++) {
        double* T_col = T.col_ptr(i);
        double tmp2 = ddot_(&d, W.col_ptr(i), &BLAS_1, T_col, &BLAS_1);
        tmp2 -= W(i, i) * T_col[i];
        double tmp1 = 1.0 / (W(i, i) - tmp2);
        double neg_tmp1 = -tmp1;
        dscal_(&d, &neg_tmp1, T_col, &BLAS_1);
        T_col[i] = tmp1;
    }
    for (int i = 0; i < d; i++) df += size_a[i];
}

// Determinant via Gaussian elimination (for log-likelihood)
static double det_gaussian(const Matrix& m) {
    const int n = m.rows;
    if (n != m.cols) return 0.0;
    double det_sign = 1.0, det_val = 1.0;
    constexpr double eps = 1e-15;
    // Row-major copy for pivot
    std::vector<double> rm(static_cast<size_t>(n) * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            rm[static_cast<size_t>(i) * n + j] = m(i, j);

    for (int k = 0; k < n; k++) {
        int pivot = k;
        double max_abs = std::fabs(rm[static_cast<size_t>(k) * n + k]);
        for (int i = k + 1; i < n; i++) {
            double av = std::fabs(rm[static_cast<size_t>(i) * n + k]);
            if (av > max_abs) { max_abs = av; pivot = i; }
        }
        if (max_abs <= eps) return 0.0;
        if (pivot != k) {
            for (int j = k; j < n; j++)
                std::swap(rm[static_cast<size_t>(k) * n + j],
                          rm[static_cast<size_t>(pivot) * n + j]);
            det_sign = -det_sign;
        }
        double diag = rm[static_cast<size_t>(k) * n + k];
        det_val *= diag;
        for (int i = k + 1; i < n; i++) {
            double factor = rm[static_cast<size_t>(i) * n + k] / diag;
            if (factor == 0.0) continue;
            for (int j = k + 1; j < n; j++)
                rm[static_cast<size_t>(i) * n + j] -= factor * rm[static_cast<size_t>(k) * n + j];
        }
    }
    return det_sign * det_val;
}

static double trace_matmul(const Matrix& a, const Matrix& b) {
    // trace(A * B) = sum_i dot(A[i,:], B[:,i]) = sum_i dot(A_col_transposed, B_col)
    // For column-major: A(i,k) = a.v[k*d+i], B(k,i) = b.v[i*d+k]
    // So trace = sum over columns: dot(a_col_k, b_col_k) summed? No.
    // trace(A*B) = sum_{i,k} A(i,k)*B(k,i) = sum_k dot(A[:,k], B[k,:])
    // = sum_k dot(col_k(A), row_k(B))
    // For column-major, col_k(A) is contiguous. row_k(B) is NOT contiguous.
    // But: trace(A*B) = sum of element-wise product of A and B^T
    // = dot(vec(A), vec(B^T)). Since B is column-major, B^T row-major = column of B^T = row of B.
    // Simpler: trace(A*B) = sum_i (A^T)[:,i] . B[:,i] = sum of A^T .* B element-wise
    // For symmetric case (our W and S are symmetric): trace = dot(vec(A), vec(B))
    // Actually even simpler for general case:
    // trace(A*B) = sum_{i} sum_{k} A(i,k)*B(k,i)
    //            = sum_{k} sum_{i} A(i,k)*B(k,i)
    //            = sum_{k} dot(A_col_k, B_row_k_as_col)
    // B_row_k = B(k,0..d-1), stride = d in column-major storage
    // A_col_k is contiguous
    const int d = a.rows;
    double tr = 0;
    for (int k = 0; k < d; k++) {
        // dot( A[:,k], B[k,:] ) where B[k,:] has stride d in memory
        tr += ddot_(&d, a.col_ptr(k), &BLAS_1, b.v.data() + k, &d);
    }
    return tr;
}

// =========================================================================
// Glasso: outer driver with pre-screening
// =========================================================================

GlassoResult glasso(const double* S_data, int d,
                    const double* lambda, int nlambda,
                    bool scr, bool cov_output)
{
    GlassoResult res;
    res.loglik.assign(nlambda, -static_cast<double>(d));
    res.sparsity.assign(nlambda, 0.0);
    res.df.assign(nlambda, 0);
    res.path.resize(nlambda);
    res.icov.resize(nlambda);
    if (cov_output) res.cov.resize(nlambda);
    if (nlambda == 0 || d <= 0) return res;

    // Build S matrix
    Matrix S(d, d);
    std::memcpy(S.v.data(), S_data, static_cast<size_t>(d) * d * sizeof(double));

    std::vector<double> s_diag(d);
    for (int i = 0; i < d; i++) s_diag[i] = S(i, i);

    const double sparsity_denom = d > 1 ? static_cast<double>(d) * (d - 1) : 1.0;
    bool zero_sol = true;

    std::vector<Matrix> tmp_icov(nlambda), tmp_path(nlambda), tmp_cov;
    if (cov_output) tmp_cov.resize(nlambda);
    const Matrix* prev_icov_ptr = nullptr;
    const Matrix* prev_cov_ptr = nullptr;
    Matrix prev_cov_buffer;

    double prev_lambda = -1.0;

    for (int i = nlambda - 1; i >= 0; i--) {
        double lambda_i = lambda[i];

        // --- Strong Rule screening for active rows ---
        // For the first lambda (no previous solution), use basic screening.
        // For subsequent lambdas, use strong rule: include row_i if any
        // off-diagonal |S(row_i,col_i) - W_prev(row_i,col_i)| >= strong_thresh,
        // where strong_thresh = 2*lambda_i - prev_lambda.
        bool use_strong_rule = (prev_cov_ptr != nullptr && prev_lambda > 0);
        double strong_thresh = use_strong_rule
            ? std::max(0.0, 2.0 * lambda_i - prev_lambda) : 0.0;

        auto compute_active_set = [&](bool strong) -> std::vector<int> {
            std::vector<int> active;
            active.reserve(d);
            for (int row_i = 0; row_i < d; row_i++) {
                bool include = false;
                if (strong) {
                    // Strong rule: check gradient gap from previous W
                    for (int col_i = 0; col_i < d; col_i++) {
                        if (col_i == row_i) continue;
                        double gap = std::fabs(S(row_i, col_i)
                                               - (*prev_cov_ptr)(row_i, col_i));
                        if (gap >= strong_thresh) { include = true; break; }
                    }
                } else {
                    // Basic screening: include if >=2 entries exceed lambda
                    int cnt = 0;
                    for (int col_i = 0; col_i < d; col_i++) {
                        if (std::fabs(S(row_i, col_i)) > lambda_i) {
                            if (++cnt > 1) { include = true; break; }
                        }
                    }
                }
                if (include) active.push_back(row_i);
            }
            return active;
        };

        std::vector<int> z = compute_active_set(use_strong_rule);

        // Solve sub-problem with KKT safety check for strong rule.
        // If strong rule wrongly excluded rows, expand z and re-solve once.
        int q = static_cast<int>(z.size());
        Matrix sub_S(q, q), sub_W(q, q), sub_T(q, q);
        int sub_df = 0;

        auto fill_and_solve = [&]() {
            q = static_cast<int>(z.size());
            sub_S.resize(q, q); sub_W.resize(q, q); sub_T.resize(q, q);
            sub_df = 0;
            for (int ii = 0; ii < q; ii++) {
                for (int jj = 0; jj < q; jj++) {
                    sub_S(ii, jj) = S(z[ii], z[jj]);
                    if (zero_sol || prev_cov_ptr == nullptr || prev_icov_ptr == nullptr) {
                        sub_W(ii, jj) = S(z[ii], z[jj]);
                        sub_T(ii, jj) = (ii == jj) ? 1.0 : 0.0;
                    } else {
                        sub_W(ii, jj) = (*prev_cov_ptr)(z[ii], z[jj]);
                        sub_T(ii, jj) = (*prev_icov_ptr)(z[ii], z[jj]);
                    }
                }
            }
            if (q > 0) {
                glasso_sub(sub_S, sub_W, sub_T, q, lambda_i, sub_df, scr);
                zero_sol = false;
            } else {
                zero_sol = true;
            }
        };

        fill_and_solve();

        // KKT check: only needed when strong rule was used
        if (use_strong_rule) {
            std::vector<bool> in_z(d, false);
            for (int idx : z) in_z[idx] = true;

            bool violation = false;
            for (int row_i = 0; row_i < d; row_i++) {
                if (in_z[row_i]) continue;
                // Excluded rows have W off-diag = 0, so KKT check
                // reduces to basic screening: |S(row_i,col_i)| > lambda_i
                int cnt = 0;
                for (int col_i = 0; col_i < d; col_i++) {
                    if (std::fabs(S(row_i, col_i)) > lambda_i) {
                        if (++cnt > 1) {
                            z.push_back(row_i);
                            violation = true;
                            break;
                        }
                    }
                }
            }
            if (violation) {
                std::sort(z.begin(), z.end());
                fill_and_solve();
            }
        }

        // Build full-size output matrices
        Matrix& cur_icov = tmp_icov[i];
        Matrix& cur_path = tmp_path[i];
        cur_icov.resize(d, d);
        cur_path.resize(d, d);

        Matrix* cur_cov_ptr;
        if (cov_output) {
            tmp_cov[i].resize(d, d);
            cur_cov_ptr = &tmp_cov[i];
        } else {
            prev_cov_buffer.resize(d, d);
            cur_cov_ptr = &prev_cov_buffer;
        }

        // Initialize diagonal
        for (int j = 0; j < d; j++) {
            cur_icov(j, j) = 1.0 / (s_diag[j] + lambda_i);
            (*cur_cov_ptr)(j, j) = s_diag[j] + lambda_i;
        }

        if (!zero_sol) {
            for (int ii = 0; ii < q; ii++) {
                for (int jj = 0; jj < q; jj++) {
                    cur_icov(z[ii], z[jj]) = sub_T(ii, jj);
                    (*cur_cov_ptr)(z[ii], z[jj]) = sub_W(ii, jj);
                    double pval = (sub_T(ii, jj) == 0) ? 0.0 : 1.0;
                    cur_path(z[ii], z[jj]) = (ii == jj) ? 0.0 : pval;
                }
            }
            res.sparsity[i] = sub_df / sparsity_denom;
            res.df[i] = sub_df / 2;
            double det_t = det_gaussian(sub_T);
            double tr = trace_matmul(sub_T, sub_S);
            res.loglik[i] = (det_t > 0) ? (std::log(det_t) - tr - (d - q))
                                         : -std::numeric_limits<double>::infinity();
        }
        prev_icov_ptr = &cur_icov;
        prev_cov_ptr = cur_cov_ptr;
        prev_lambda = lambda_i;
    }

    res.path = std::move(tmp_path);
    res.icov = std::move(tmp_icov);
    if (cov_output) res.cov = std::move(tmp_cov);
    return res;
}

// =========================================================================
// MB graph estimation (lossless screening)
// =========================================================================

MBResult mb(const double* S_data, int d,
            const double* lambda, int nlambda)
{
    MBResult res;
    if (d <= 0 || nlambda <= 0) { res.columns.resize(d > 0 ? d : 0); return res; }
    res.columns.resize(d);

    const double thol = 1e-4;
    const int MAX_ITER = 10000;

    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
    std::vector<double> w0(d, 0.0), w1(d, 0.0);
    std::vector<int> idx_a(d), idx_i(d);

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for (int m = 0; m < d; m++) {
        ColResult& col = res.columns[m];
        idx_i[m] = 0;
        for (int j = 0; j < m; j++) idx_i[j] = 1;
        for (int j = m + 1; j < d; j++) idx_i[j] = 1;

        int size_a = 0;
        std::fill(w0.begin(), w0.end(), 0.0);
        std::fill(w1.begin(), w1.end(), 0.0);

        for (int i = 0; i < nlambda; i++) {
            double ilambda = lambda[i];
            int gap_ext = 1, iter_ext = 0;
            while (gap_ext != 0 && iter_ext < MAX_ITER) {
                int size_a_prev = size_a;
                for (int j = 0; j < d; j++) {
                    if (idx_i[j] == 1) {
                        double r = cm(S_data, d, m, j);
                        for (int k = 0; k < size_a; k++)
                            r -= cm(S_data, d, j, idx_a[k]) * w0[idx_a[k]];

                        w1[j] = threshold_l1(r, ilambda);
                        if (w1[j] != 0) {
                            idx_a[size_a++] = j;
                            idx_i[j] = 0;
                        }
                        w0[j] = w1[j];
                    }
                }
                gap_ext = size_a - size_a_prev;

                double gap_int = 1;
                int iter_int = 0;
                while (gap_int > thol && iter_int < MAX_ITER) {
                    double tmp1 = 0, tmp2 = 0;
                    for (int j = 0; j < size_a; j++) {
                        int w_idx = idx_a[j];
                        double r = cm(S_data, d, m, w_idx) + w0[w_idx];
                        for (int k = 0; k < size_a; k++)
                            r -= cm(S_data, d, w_idx, idx_a[k]) * w0[idx_a[k]];

                        w1[w_idx] = threshold_l1(r, ilambda);
                        tmp2 += std::fabs(w1[w_idx]);
                        tmp1 += std::fabs(w1[w_idx] - w0[w_idx]);
                        w0[w_idx] = w1[w_idx];
                    }
                    gap_int = (tmp2 > 0) ? tmp1 / tmp2 : 0;
                    iter_int++;
                }
                int junk_a = 0;
                for (int j = 0; j < size_a; j++) {
                    int w_idx = idx_a[j];
                    if (w1[w_idx] == 0) { junk_a++; idx_i[w_idx] = 1; }
                    else idx_a[j - junk_a] = w_idx;
                }
                size_a -= junk_a;
                iter_ext++;
            }
            for (int j = 0; j < size_a; j++) {
                col.vals.push_back(w1[idx_a[j]]);
                col.indices.push_back(i * d + idx_a[j]);
            }
        }
    }
    #ifdef _OPENMP
    }
    #endif
    return res;
}

// =========================================================================
// MB graph estimation (lossy screening)
// =========================================================================

MBResult mb_scr(const double* S_data, int d,
                const double* lambda, int nlambda,
                const int* idx_scr, int nscr)
{
    MBResult res;
    if (d <= 0 || nlambda <= 0 || nscr < 0) { res.columns.resize(d > 0 ? d : 0); return res; }
    res.columns.resize(d);

    const double thol = 1e-4;
    const int MAX_ITER = 10000;

    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
    std::vector<double> w0(d, 0.0), w1(d, 0.0);
    std::vector<int> idx_a(nscr), idx_i_local(nscr);

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for (int m = 0; m < d; m++) {
        ColResult& col = res.columns[m];
        int size_a = 0;

        for (int j = 0; j < nscr; j++)
            idx_i_local[j] = idx_scr[static_cast<size_t>(m) * nscr + j];

        std::fill(w0.begin(), w0.end(), 0.0);
        std::fill(w1.begin(), w1.end(), 0.0);

        for (int i = 0; i < nlambda; i++) {
            double ilambda = lambda[i];
            int gap_ext = 1, iter_ext = 0;
            while (iter_ext < MAX_ITER && gap_ext > 0) {
                int size_a_prev = size_a;
                for (int j = 0; j < nscr; j++) {
                    int w_idx = idx_i_local[j];
                    if (w_idx != -1) {
                        double r = cm(S_data, d, m, w_idx);
                        for (int k = 0; k < size_a; k++)
                            r -= cm(S_data, d, w_idx, idx_a[k]) * w0[idx_a[k]];

                        w1[w_idx] = threshold_l1(r, ilambda);
                        if (w1[w_idx] != 0) {
                            idx_a[size_a++] = w_idx;
                            idx_i_local[j] = -1;
                        }
                        w0[w_idx] = w1[w_idx];
                    }
                }
                gap_ext = size_a - size_a_prev;

                double gap_int = 1;
                int iter_int = 0;
                while (gap_int > thol && iter_int < MAX_ITER) {
                    double tmp1 = 0, tmp2 = 0;
                    for (int j = 0; j < size_a; j++) {
                        int w_idx = idx_a[j];
                        double r = cm(S_data, d, m, w_idx) + w0[w_idx];
                        for (int k = 0; k < size_a; k++)
                            r -= cm(S_data, d, w_idx, idx_a[k]) * w0[idx_a[k]];

                        w1[w_idx] = threshold_l1(r, ilambda);
                        tmp2 += std::fabs(w1[w_idx]);
                        tmp1 += std::fabs(w1[w_idx] - w0[w_idx]);
                        w0[w_idx] = w1[w_idx];
                    }
                    gap_int = (tmp2 > 0) ? tmp1 / tmp2 : 0;
                    iter_int++;
                }
                iter_ext++;
            }
            for (int j = 0; j < size_a; j++) {
                col.vals.push_back(w1[idx_a[j]]);
                col.indices.push_back(i * d + idx_a[j]);
            }
        }
    }
    #ifdef _OPENMP
    }
    #endif
    return res;
}

// =========================================================================
// TIGER (sqrt-lasso graph estimation)
// =========================================================================

TigerResult tiger(const double* data_colmajor, int n, int d,
                  const double* lambda, int nlambda)
{
    TigerResult res;
    if (n <= 0 || d <= 0 || nlambda <= 0) {
        res.columns.resize(d > 0 ? d : 0);
        return res;
    }
    res.columns.resize(d);

    // Pre-initialize icov matrices
    res.icov.resize(nlambda);
    for (int i = 0; i < nlambda; i++) res.icov[i].resize(d, d);

    const double prec = 1e-4;
    const int max_iter = 1000;
    const int num_relaxation_round = 3;
    const double eps = 1e-12;

    // Accessor for data (column-major)
    auto x_at = [&](int r, int c) -> double {
        return data_colmajor[static_cast<size_t>(c) * n + r];
    };

    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
    std::vector<double> Xb(n, 0.0), r_vec(n, 0.0), grad(d, 0.0), w1(d, 0.0);
    std::vector<double> Y(n, 0.0), gr(d, 0.0);
    std::vector<double> Xb_master(n, 0.0), w1_master(d, 0.0);
    std::vector<int> actset_indcat(d, 0), actset_indcat_master(d, 0);
    std::vector<int> actset_idx;
    std::vector<double> old_coef(d, 0.0), grad_master(d, 0.0);

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for (int m = 0; m < d; m++) {
        ColResult& col = res.columns[m];

        std::fill(Xb.begin(), Xb.end(), 0.0);
        std::fill(w1.begin(), w1.end(), 0.0);
        std::fill(grad.begin(), grad.end(), 0.0);
        std::fill(gr.begin(), gr.end(), 0.0);
        std::fill(w1_master.begin(), w1_master.end(), 0.0);
        std::fill(Xb_master.begin(), Xb_master.end(), 0.0);
        std::fill(actset_indcat.begin(), actset_indcat.end(), 0);
        std::fill(actset_indcat_master.begin(), actset_indcat_master.end(), 0);
        std::fill(old_coef.begin(), old_coef.end(), 0.0);
        std::fill(grad_master.begin(), grad_master.end(), 0.0);

        for (int t = 0; t < n; t++) Y[t] = x_at(t, m);

        double L = 0, sum_r2 = 0;
        auto refresh_residual = [&]() {
            sum_r2 = 0;
            for (int t = 0; t < n; t++) sum_r2 += r_vec[t] * r_vec[t];
            if (sum_r2 < eps) sum_r2 = eps;
            L = std::sqrt(sum_r2 / n);
            if (L < eps) L = eps;
        };

        for (int t = 0; t < n; t++) r_vec[t] = Y[t] - Xb[t];
        refresh_residual();
        double dev_thr = std::fabs(L) * prec;

        for (int j = 0; j < d; j++) {
            const double* x_col = data_colmajor + static_cast<size_t>(j) * n;
            double s = ddot_(&n, r_vec.data(), &BLAS_1, x_col, &BLAS_1);
            grad[j] = s / (n * L);
            gr[j] = std::fabs(grad[j]);
            grad_master[j] = gr[j];
            w1_master[j] = w1[j];
        }
        std::memcpy(Xb_master.data(), Xb.data(), n * sizeof(double));

        for (int i = 0; i < nlambda; i++) {
            double stage_lambda = lambda[i];
            w1 = w1_master;
            Xb = Xb_master;
            for (int j = 0; j < d; j++) {
                gr[j] = grad_master[j];
                actset_indcat[j] = actset_indcat_master[j];
            }

            double threshold = (i > 0) ? 2 * lambda[i] - lambda[i - 1] : 2 * lambda[i];
            for (int j = 0; j < d; j++)
                if (j != m && gr[j] > threshold) actset_indcat[j] = 1;

            for (int t = 0; t < n; t++) r_vec[t] = Y[t] - Xb[t];
            refresh_residual();

            double tmp_change = 0, local_change = 0;

            auto update_coordinate = [&](int ci) {
                const double* x_col = data_colmajor + static_cast<size_t>(ci) * n;
                double sum_wxx = 0, sum_g = 0;
                double inv_sum_r2 = 1.0 / sum_r2;
                for (int t = 0; t < n; t++) {
                    double xv = x_col[t];
                    double rv = r_vec[t];
                    double wxx_val = (1.0 - rv * rv * inv_sum_r2) * xv * xv;
                    sum_wxx += wxx_val;
                    sum_g += wxx_val * w1[ci] + rv * xv;
                }
                double g = sum_g / (n * L);
                double a = sum_wxx / (n * L);
                double oldv = w1[ci];
                w1[ci] = (std::fabs(a) > eps) ? threshold_l1(g, stage_lambda) / a : 0.0;
                double delta = w1[ci] - oldv;
                if (delta != 0) {
                    // Xb += delta * X[:,ci];  r_vec -= delta * X[:,ci]
                    daxpy_(&n, &delta, x_col, &BLAS_1, Xb.data(), &BLAS_1);
                    double neg_delta = -delta;
                    daxpy_(&n, &neg_delta, x_col, &BLAS_1, r_vec.data(), &BLAS_1);
                    refresh_residual();
                }
            };

            int loopcnt_level_0 = 0;
            while (loopcnt_level_0 < num_relaxation_round) {
                loopcnt_level_0++;

                int loopcnt_level_1 = 0;
                while (loopcnt_level_1 < max_iter) {
                    loopcnt_level_1++;
                    for (int j = 0; j < d; j++) old_coef[j] = w1[j];
                    refresh_residual();

                    actset_idx.clear();
                    for (int j = 0; j < d; j++) {
                        if (j == m || !actset_indcat[j]) continue;
                        update_coordinate(j);
                        if (std::fabs(w1[j]) > 0) actset_idx.push_back(j);
                    }

                    // Level 2: proximal newton on active set
                    int loopcnt_level_2 = 0;
                    while (loopcnt_level_2 < max_iter) {
                        loopcnt_level_2++;
                        bool term2 = true;
                        for (size_t k = 0; k < actset_idx.size(); k++) {
                            int idx = actset_idx[k];
                            double old_w1 = w1[idx];
                            update_coordinate(idx);
                            tmp_change = old_w1 - w1[idx];
                            const double* xc = data_colmajor + static_cast<size_t>(idx) * n;
                            double hsum = 0;
                            double inv_LLn = 1.0 / (L * L * n);
                            for (int t = 0; t < n; t++) {
                                double xv = xc[t];
                                double rv = r_vec[t];
                                hsum += xv * xv * (1.0 - rv * rv * inv_LLn);
                            }
                            double h = std::fabs(hsum / (n * L));
                            local_change = h * tmp_change * tmp_change / (2.0 * L * n);
                            if (local_change > dev_thr) term2 = false;
                        }
                        if (term2) break;
                    }

                    // Check stopping criterion 1
                    bool term1 = true;
                    for (size_t k = 0; k < actset_idx.size(); k++) {
                        int idx = actset_idx[k];
                        tmp_change = old_coef[idx] - w1[idx];
                        const double* xc = data_colmajor + static_cast<size_t>(idx) * n;
                        double hsum = 0;
                        double inv_LLn = 1.0 / (L * L * n);
                        for (int t = 0; t < n; t++) {
                            double xv = xc[t];
                            double rv = r_vec[t];
                            hsum += xv * xv * (1.0 - rv * rv * inv_LLn);
                        }
                        double h = std::fabs(hsum / (n * L));
                        local_change = h * tmp_change * tmp_change / (2.0 * L * n);
                        if (local_change > dev_thr) term1 = false;
                    }
                    for (int t = 0; t < n; t++) r_vec[t] = Y[t] - Xb[t];
                    refresh_residual();
                    if (term1) break;

                    // Check stopping criterion 2: active set change
                    bool new_active = false;
                    for (int k = 0; k < d; k++) {
                        if (k == m || actset_indcat[k] != 0) continue;
                        const double* x_col_k = data_colmajor + static_cast<size_t>(k) * n;
                        double s = ddot_(&n, r_vec.data(), &BLAS_1, x_col_k, &BLAS_1);
                        grad[k] = s / (n * L);
                        gr[k] = std::fabs(grad[k]);
                        if (gr[k] > stage_lambda) {
                            actset_indcat[k] = 1;
                            new_active = true;
                        }
                    }
                    if (!new_active) break;
                }

                if (loopcnt_level_0 == 1) {
                    w1_master = w1;
                    Xb_master = Xb;
                    for (int j = 0; j < d; j++) {
                        grad_master[j] = gr[j];
                        actset_indcat_master[j] = actset_indcat[j];
                    }
                }
            }

            // Collect results
            for (size_t j = 0; j < actset_idx.size(); j++) {
                int w_idx = actset_idx[j];
                col.vals.push_back(w1[w_idx]);
                col.indices.push_back(i * d + w_idx);
            }

            for (int t = 0; t < n; t++) r_vec[t] = Y[t] - Xb[t];
            refresh_residual();
            double tal = L;

            // Write icov column m — each thread writes different column, no race
            Matrix& icov_ref = res.icov[i];
            icov_ref(m, m) = (tal > 0) ? 1.0 / (tal * tal) : 0.0;
            for (int j = 0; j < d; j++)
                if (j != m) icov_ref(j, m) = -icov_ref(m, m) * w1[j];
        }
    }
    #ifdef _OPENMP
    }
    #endif

    // Symmetrize icov
    for (int i = 0; i < nlambda; i++) {
        Matrix& ic = res.icov[i];
        for (int c0 = 0; c0 < d; c0++)
            for (int r0 = 0; r0 < d; r0++)
                ic(r0, c0) = 0.5 * (ic(r0, c0) + ic(c0, r0));
    }
    return res;
}

// =========================================================================
// RIC (Rotation Information Criterion)
// =========================================================================

double ric(const double* X_data, int n, int d, const int* r, int t)
{
    if (d <= 1 || n <= 0 || t <= 0) return 0.0;

    double lambda_min = std::numeric_limits<double>::infinity();

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) reduction(min:lambda_min)
    #endif
    for (int i = 0; i < t; i++) {
        int tmp_r = r[i];
        if (tmp_r < 0) tmp_r = 0;
        if (tmp_r > n) tmp_r = n;
        int split = n - tmp_r;
        double lambda_max = 0;
        for (int j = 0; j < d; j++) {
            for (int k = j + 1; k < d; k++) {
                double tmp = 0;
                for (int mm = 0; mm < split; mm++)
                    tmp += cm(X_data, n, mm + tmp_r, j) * cm(X_data, n, mm, k);
                for (int mm = split; mm < n; mm++)
                    tmp += cm(X_data, n, mm - split, j) * cm(X_data, n, mm, k);
                tmp = std::fabs(tmp);
                if (tmp > lambda_max) lambda_max = tmp;
            }
        }
        if (lambda_max < lambda_min) lambda_min = lambda_max;
    }
    if (!std::isfinite(lambda_min)) return 0.0;
    return lambda_min;
}

// =========================================================================
// Scale-free graph generator
// =========================================================================

void sfgen(int d0, int d, int* G_out, const double* rands)
{
    // G_out: d*d column-major, pre-zeroed by caller
    std::memset(G_out, 0, static_cast<size_t>(d) * d * sizeof(int));
    std::vector<int> degree(d, 0);

    // Initial cycle of d0 nodes
    for (int i = 0; i < d0 - 1; i++) {
        G_out[static_cast<size_t>(i) * d + (i + 1)] = 1;
        G_out[static_cast<size_t>(i + 1) * d + i] = 1;
    }
    G_out[static_cast<size_t>(0) * d + (d0 - 1)] = 1;
    G_out[static_cast<size_t>(d0 - 1) * d + 0] = 1;

    for (int i = 0; i < d0; i++) degree[i] = 2;
    int total = 2 * d0;

    for (int i = d0; i < d; i++) {
        double x = static_cast<double>(total) * rands[i - d0];
        int tmp = 0, j = 0;
        while (tmp < x && j < i) { tmp += degree[j]; j++; }
        if (j > 0) j--;
        G_out[static_cast<size_t>(i) * d + j] = 1;
        G_out[static_cast<size_t>(j) * d + i] = 1;
        total += 2;
        degree[j]++;
        degree[i]++;
    }
}

} // namespace huge
