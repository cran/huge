#include "math.h"
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

void hugeglasso_sub(Eigen::MatrixXd &S, Eigen::MatrixXd &W, Eigen::MatrixXd &T, int d, double ilambda, int &df, bool scr);

//[[Rcpp::export]]
List hugeglasso(Eigen::Map<Eigen::MatrixXd> S, NumericVector lambda, bool scr, bool verbose, bool cov_output)
{
    unsigned int nlambda = lambda.size();
    int nfeatures = S.rows();
    auto S_diag = S.diagonal().array();
    int &d = nfeatures;
    double sparsity_denom = d > 1 ? static_cast<double>(d) * (d - 1) : 1.0;

    bool zero_sol = true;

    // define result list: loglik
    NumericVector loglik(nlambda,-nfeatures), sparsity(nlambda);
    IntegerVector df(nlambda);
    List result;
    result["loglik"] = loglik;
    result["sparsity"] = sparsity;
    result["df"] = df;
    if(nlambda == 0){
        result["path"] = List::create();
        result["icov"] = List::create();
        if(cov_output) result["cov"] = List::create();
        return result;
    }

    std::vector<Eigen::MatrixXd> tmp_icov_p(nlambda), tmp_path_p(nlambda), tmp_cov_p;
    if(cov_output) tmp_cov_p.resize(nlambda);
    const Eigen::MatrixXd *prev_icov_ptr = nullptr;
    const Eigen::MatrixXd *prev_cov_ptr = nullptr;
    Eigen::MatrixXd prev_cov_buffer;
    for(int i = nlambda-1; i >= 0; i--)
    {
        double lambda_i = lambda[i];
        // pre-screening by z
        vector<int> z;
        z.reserve(d);
        for (int row_i = 0; row_i < d; row_i++) {
            int break_flag = 0;
            for (int col_i = 0; col_i < d; col_i++) {
                if(break_flag > 1) break;
                if(S(row_i,col_i) > lambda_i or S(row_i,col_i) < -lambda_i) break_flag++;
            }
            if(break_flag > 1) z.push_back(row_i);
        }
        int q = z.size();
        MatrixXd sub_S(q,q), sub_W(q,q), sub_T(q,q);
        int sub_df=0;

        //#pragma omp parallel for
        for (int ii = 0; ii < q; ii++) {
            for (int jj = 0; jj < q; jj++) {
                sub_S(ii,jj) = S(z[ii],z[jj]);
                if(zero_sol || prev_cov_ptr == nullptr || prev_icov_ptr == nullptr) {
                    sub_W(ii,jj) = S(z[ii],z[jj]);
                    sub_T(ii,jj) = ii==jj ? 1 : 0;
                }
                else {
                    sub_W(ii,jj) = (*prev_cov_ptr)(z[ii],z[jj]);
                    sub_T(ii,jj) = (*prev_icov_ptr)(z[ii],z[jj]);
                }
            }
        }

        if(q>0)
        {
            if(verbose){
              if(scr)
                Rcout << "\rConducting the graphical lasso (glasso) wtih lossy screening....in progress: " << floor(100*(1-1.*i/nlambda))<<"%";
              if(!scr)
                Rcout << "\rConducting the graphical lasso (glasso) wtih lossless screening....in progress: " << floor(100*(1-1.*i/nlambda))<<"%";
            }

            hugeglasso_sub(sub_S, sub_W, sub_T, q, lambda_i, sub_df, scr);
            zero_sol = false;
        }
        if(q == 0) zero_sol = true;
        // update result list
        Eigen::MatrixXd &tmp_icov = tmp_icov_p[i];
        Eigen::MatrixXd &tmp_path = tmp_path_p[i];
        tmp_icov.resize(d, d);
        tmp_path.resize(d, d);
        Eigen::MatrixXd *tmp_cov_ptr;
        if(cov_output){
            Eigen::MatrixXd &tmp_cov_ref = tmp_cov_p[i];
            tmp_cov_ref.resize(d, d);
            tmp_cov_ptr = &tmp_cov_ref;
        } else {
            prev_cov_buffer.resize(d, d);
            tmp_cov_ptr = &prev_cov_buffer;
        }
        tmp_icov.setZero();
        tmp_icov.diagonal() = 1/(S_diag + lambda_i);
        tmp_cov_ptr->setZero();
        tmp_cov_ptr->diagonal() = S_diag + lambda_i;
        tmp_path.setZero();

        if(!zero_sol)
        {
            //#pragma omp parallel for
            for (int ii = 0; ii < q; ii++) {
                for (int jj = 0; jj < q; jj++) {
                    tmp_icov(z[ii],z[jj]) = sub_T(ii,jj);
                    (*tmp_cov_ptr)(z[ii],z[jj]) = sub_W(ii,jj);
                    tmp_path(z[ii],z[jj]) = sub_T(ii,jj)==0 ? 0 : 1;
                    tmp_path(z[ii],z[jj]) = ii==jj ? 0 : tmp_path(z[ii],z[jj]);
                }
            }
            sparsity[i] = sub_df / sparsity_denom;
            df[i] = sub_df/2;
            loglik[i] = log(sub_T.determinant()) - (sub_T * sub_S).diagonal().sum() - (d-q);
        }
        prev_icov_ptr = &tmp_icov;
        prev_cov_ptr = tmp_cov_ptr;
    }

    List path(nlambda), icov(nlambda), cov;
    if(cov_output) cov = List(nlambda);
    for(unsigned int i = 0; i < nlambda; i++){
        path[i] = tmp_path_p[i];
        icov[i] = tmp_icov_p[i];
        if(cov_output) cov[i] = tmp_cov_p[i];
    }
    result["path"] = path;
    result["icov"] = icov;
    if(cov_output) result["cov"] = cov;
    return result;
}

void hugeglasso_sub(Eigen::MatrixXd &S, Eigen::MatrixXd &W, Eigen::MatrixXd &T, int d, double ilambda, int &df, bool scr)
{
    int i,j,k;
    int rss_idx,w_idx;

    int gap_int;
    double gap_ext,gap_act;
    const double thol_act = 1e-4;
    const double thol_ext = 1e-4;

    const int MAX_ITER_EXT = 100;
    const int MAX_ITER_INT = 10000;
    const int MAX_ITER_ACT = 10000;
    int iter_ext,iter_int,iter_act;

    Eigen::MatrixXi idx_a(d, d); // active set
    Eigen::MatrixXi idx_i(d, d); // The set possibly can join active set
    std::vector<int> size_a(d, 0); //sizes of active sets
    std::vector<double> w1(d, 0.0);
    std::vector<double> ww(d, 0.0);

    int size_a_prev;
    int junk_a;

    double r;
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    const double neg_ilambda = -ilambda;

    // recover initial solution for each individual lasso
    for(i=0;i<d;i++){
        int *idx_a_col = idx_a.col(i).data();
        int *idx_i_col = idx_i.col(i).data();
        double *T_col = T.col(i).data();
        const double *S_col = S.col(i).data();

        W(i, i) = S_col[i] + ilambda;
        size_a[i] = 0;
        tmp1 = T_col[i];
        T_col[i] = 0;

        for(j=0;j<d;j++){
            if(scr && fabs(S_col[j]) <= ilambda){
                idx_i_col[j] = -1;
                T_col[j] = 0;
                continue;
            }

            if(T_col[j] != 0){
                idx_a_col[size_a[i]++] = j;
                idx_i_col[j] = -1;
                T_col[j] = -T_col[j]/tmp1;
            } else {
                idx_i_col[j] = 1;
            }
        }
        idx_i_col[i] = -1;
    }

    gap_ext = 1;
    iter_ext = 0;
    while(gap_ext > thol_ext && iter_ext < MAX_ITER_EXT)
    {
        tmp6 = 0;
        tmp5 = 0;
        for(i=0;i<d;i++)
        {
            int active_size = size_a[i];
            int *idx_a_col = idx_a.col(i).data();
            int *idx_i_col = idx_i.col(i).data();
            double *T_col = T.col(i).data();
            const double *S_col = S.col(i).data();

            gap_int = 1;
            iter_int = 0;

            for(j=0;j<d;j++)
                ww[j] = T_col[j];

            while(gap_int != 0 && iter_int < MAX_ITER_INT)
            {
                size_a_prev = active_size;
                for(j=0;j<d;j++)
                {
                    if(idx_i_col[j] != -1)
                    {
                        r = S_col[j];
                        for(k=0;k<active_size;k++)
                        {
                            rss_idx = idx_a_col[k];
                            r -= W(rss_idx, j) * T_col[rss_idx];
                        }

                        double w_new = 0.0;
                        if(r > ilambda)
                        {
                            w_new = (r - ilambda)/W(j, j);
                            idx_a_col[active_size++] = j;
                            idx_i_col[j] = -1;
                        }
                        else if(r < neg_ilambda)
                        {
                            w_new = (r + ilambda)/W(j, j);
                            idx_a_col[active_size++] = j;
                            idx_i_col[j] = -1;
                        }

                        w1[j] = w_new;
                        T_col[j] = w_new;
                    }
                }
                gap_int = active_size - size_a_prev;

                gap_act = 1;
                iter_act = 0;

                while(gap_act > thol_act && iter_act < MAX_ITER_ACT)
                {
                    tmp3 = 0;
                    tmp4 = 0;
                    for(j=0;j<active_size;j++)
                    {
                        w_idx = idx_a_col[j];
                        r = S_col[w_idx] + T_col[w_idx] * W(w_idx, w_idx);
                        for(k=0;k<active_size;k++)
                        {
                            rss_idx = idx_a_col[k];
                            r -= W(rss_idx, w_idx) * T_col[rss_idx];
                        }

                        double w_new = 0.0;
                        if(r > ilambda)
                            w_new = (r - ilambda)/W(w_idx, w_idx);
                        else if(r < neg_ilambda)
                            w_new = (r + ilambda)/W(w_idx, w_idx);

                        tmp4 += fabs(w_new);
                        tmp3 += fabs(w_new - T_col[w_idx]);
                        w1[w_idx] = w_new;
                        T_col[w_idx] = w_new;
                    }
                    gap_act = tmp4 > 0 ? tmp3/tmp4 : 0;
                    iter_act++;
                }

                // move false active variables back to inactive set
                junk_a = 0;
                for(j=0;j<active_size;j++){
                    w_idx = idx_a_col[j];
                    if(w1[w_idx] == 0){
                        junk_a++;
                        idx_i_col[w_idx] = 1;
                    } else {
                        idx_a_col[j-junk_a] = w_idx;
                    }
                }
                active_size -= junk_a;
                iter_int++;
            }
            size_a[i] = active_size;

            // update W from current T column
            Eigen::Map<Eigen::VectorXd> T_col_map(T_col, d);
            Eigen::VectorXd temp = W * T_col_map;
            for(j=0;j<i;j++){
                W(j, i) = temp(j);
                W(i, j) = temp(j);
            }
            for(j=i+1;j<d;j++){
                W(j, i) = temp(j);
                W(i, j) = temp(j);
            }

            for(j=0;j<d;j++)
                tmp5 += fabs(ww[j] - T_col[j]);
            tmp6 += tmp4;
        }
        gap_ext = tmp6 > 0 ? tmp5/tmp6 : 0;
        iter_ext++;
    }

    for(i=0;i<d;i++) //Compute the final T
    {
        double *T_col = T.col(i).data();
        Eigen::Map<Eigen::VectorXd> T_col_map(T_col, d);
        tmp2 = W.col(i).dot(T_col_map) - W(i, i)*T_col[i];

        tmp1 = 1/(W(i, i)-tmp2);
        T_col_map *= -tmp1;
        T_col[i] = tmp1;
    }
    for(i=0;i<d;i++)
        df += size_a[i];
}
