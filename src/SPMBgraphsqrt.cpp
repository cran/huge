#include "math.h"
#include <vector>
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]


double thresholdl1(double x, double thr);

//[[Rcpp::export]]
List SPMBgraphsqrt(Eigen::Map<Eigen::MatrixXd> data, NumericVector lambda, int nlambda, int d)
{

    Eigen::ArrayXd Xb, r, grad, w1, Y, gr;
    Eigen::ArrayXXd X;
    Eigen::MatrixXd tmp_icov;
    tmp_icov.resize(d, d);
    tmp_icov.setZero();
    std::vector<Eigen::MatrixXd > tmp_icov_p;
    tmp_icov_p.clear();
    tmp_icov_p.reserve(nlambda > 0 ? nlambda : 0);
    for(int i = 0; i < nlambda; i++)
      tmp_icov_p.push_back(tmp_icov);
    int n = data.rows();
    if(n <= 0 || d <= 0 || nlambda <= 0)
    {
      int d_safe = d > 0 ? d : 0;
      return List::create(
        _["col_cnz"] = IntegerVector(d_safe + 1),
        _["row_idx"] = IntegerVector(0),
        _["x"] = NumericVector(0),
        _["icov"] = std::vector<Eigen::MatrixXd>()
      );
    }
    int maxdf = 0;
    maxdf = (n < d ? n : d);
    NumericVector x(d*maxdf*nlambda);
    IntegerVector col_cnz(d+1);
    IntegerVector row_idx(d*maxdf*nlambda);
    X = data;
    double prec = 1e-4;
    int max_iter = 1000;
    int num_relaxation_round = 3;
	  int cnz = 0;
    for(int m=0;m<d;m++)
    {
      Xb.resize(n);
      Xb.setZero();
      grad.resize(d);
      grad.setZero();
      gr.resize(d);
      gr.setZero();
      w1.resize(d);
      w1.setZero();
      r.resize(n);
      r.setZero();
      Y = X.col(m);

      Eigen::ArrayXd Xb_master(n);
      Eigen::ArrayXd w1_master(d);
      std::vector<int> actset_indcat(d, 0);
      std::vector<int> actset_indcat_master(d, 0);
      std::vector<int> actset_idx;
      std::vector<double> old_coef(d, 0);
      std::vector<double> grad_master(d, 0);

      double L = 0, sum_r2 = 0;
      double tmp_change = 0, local_change = 0;
      const double eps = 1e-12;
      auto refresh_residual = [&]() {
        sum_r2 = r.matrix().dot(r.matrix());
        if (sum_r2 < eps)
          sum_r2 = eps;
        L = sqrt(sum_r2 / n);
        if (L < eps)
          L = eps;
      };

      r = Y - Xb;
      refresh_residual();

      double dev_thr = fabs(L) * prec;

      //cout<<dev_thr<<endl;


      for(int i = 0; i < d; i++)
      {
        grad[i] = (r * X.col(i)).sum() / (n*L);
      }
      for(int i = 0; i < d; i++) gr[i] = fabs(grad[i]);
      w1_master = w1;
      Xb_master = Xb;
      for (int i = 0; i < d; i++) grad_master[i] = gr[i];

      for(int i=0;i<nlambda;i++)
      {
        const double stage_lambda = lambda[i];
        w1 = w1_master;
        Xb = Xb_master;

        for (int j = 0; j < d; j++)
        {
          gr[j] = grad_master[j];
          actset_indcat[j] = actset_indcat_master[j];
        }

        // init the active set
        double threshold;
        if (i > 0)
          threshold = 2 * lambda[i] - lambda[i - 1];
        else
          threshold = 2 * lambda[i];

        for (int j = 0; j < d; ++j)
          if (j != m && gr[j] > threshold)
            actset_indcat[j] = 1;
        r = Y - Xb;
        refresh_residual();
        // loop level 0: multistage convex relaxation
        int loopcnt_level_0 = 0;
        int idx;
        double old_w1;
        while (loopcnt_level_0 < num_relaxation_round)
        {
          loopcnt_level_0++;

          // loop level 1: active set update
          int loopcnt_level_1 = 0;
          bool terminate_loop_level_1 = true;
          Eigen::ArrayXd wXX(n);
          while (loopcnt_level_1 < max_iter)
          {
            loopcnt_level_1++;
            terminate_loop_level_1 = true;

            for (int j = 0; j < d; j++) old_coef[j] = w1[j];
            refresh_residual();

            // initialize actset_idx
            actset_idx.clear();
            auto update_coordinate = [&](int coord_idx) {
              wXX = (1 - r*r/sum_r2) * X.col(coord_idx) * X.col(coord_idx);
              double g = (wXX * w1[coord_idx] + r * X.col(coord_idx)).sum()/(n*L);
              double a = wXX.sum()/(n*L);
              double tmp = w1[coord_idx];
              if (fabs(a) > eps)
                w1[coord_idx] = thresholdl1(g, stage_lambda) / a;
              else
                w1[coord_idx] = 0;

              tmp = w1[coord_idx] - tmp;
              if (tmp != 0) {
                Xb = Xb + tmp * X.col(coord_idx);
                r = r - tmp * X.col(coord_idx);
                refresh_residual();
              }
            };

            for (int j = 0; j < d; ++j) {
              if (j == m || !actset_indcat[j])
                continue;
              update_coordinate(j);
              if (fabs(w1[j]) > 0)
                actset_idx.push_back(j);
            }

              // loop level 2: proximal newton on active set
            int loopcnt_level_2 = 0;
            bool terminate_loop_level_2 = true;
            while (loopcnt_level_2 < max_iter)
            {
              loopcnt_level_2++;
              terminate_loop_level_2 = true;

              for (unsigned int k = 0; k < actset_idx.size(); k++)
              {
                  idx = actset_idx[k];

                  old_w1 = w1[idx];
                  update_coordinate(idx);
                  tmp_change = old_w1 - w1[idx];
                  double h = fabs((X.col(idx) * X.col(idx) * (1 - r * r/(L*L*n))).sum()/(n*L));
                  local_change = h * tmp_change * tmp_change / (2 * L * n);
                  if (local_change > dev_thr)
                    terminate_loop_level_2 = false;
                }
              if (terminate_loop_level_2)
                break;
            }

            terminate_loop_level_1 = true;
              // check stopping criterion 1: fvalue change
            for (unsigned int k = 0; k < actset_idx.size(); ++k)
            {
              idx = actset_idx[k];
              tmp_change = old_coef[idx] - w1[idx];
              double h = fabs((X.col(idx) * X.col(idx) * (1 - r * r/(L*L*n))).sum()/(n*L));
              local_change = h * tmp_change * tmp_change / (2 * L * n);
              if (local_change > dev_thr)
                terminate_loop_level_1 = false;
            }

            r = Y - Xb;
            refresh_residual();

            if (terminate_loop_level_1)
              break;


              // check stopping criterion 2: active set change
            bool new_active_idx = false;
            for (int k = 0; k < d; ++k)
              if (k != m && actset_indcat[k] == 0)
              {
                grad[k] = (r * X.col(k)).sum() / (n*L);
                gr[k] = fabs(grad[k]);
                if (gr[k] > stage_lambda)
                {
                  actset_indcat[k] = 1;
                  new_active_idx = true;
                }
              }
            if(!new_active_idx)
              break;
          }
          if (loopcnt_level_0 == 1)
          {
            for (int j = 0; j < d; j++)
            {
              w1_master[j] = w1[j];

              grad_master[j] = gr[j];
              actset_indcat_master[j] = actset_indcat[j];
            }

            for (int j = 0; j < n; j++) Xb_master[j] = Xb[j];
          }
        }
        for(unsigned int j=0;j<actset_idx.size();j++)
        {
          int w_idx = actset_idx[j];
          x[cnz] = w1[w_idx];
          row_idx[cnz] = i*d+w_idx;
          cnz++;
          //cout<<cnz<<"    ";
        }
        r = Y - Xb;
        refresh_residual();
        double tal = L;

        Eigen::MatrixXd &tmp_icov_ref = tmp_icov_p[i];
        if(tal > 0)
          tmp_icov_ref(m, m) = 1.0/(tal*tal);
        else
          tmp_icov_ref(m, m) = 0;
        for(int j = 0; j < d; ++j)
          if (j != m)
            tmp_icov_ref(j, m) = -tmp_icov_ref(m, m)*w1[j];
      }
      col_cnz[m+1]=cnz;
    }
    for(int i = 0; i < nlambda; i++)
      tmp_icov_p[i] = (tmp_icov_p[i].transpose()+tmp_icov_p[i])/2;
	  return List::create(
	    _["col_cnz"] = col_cnz,
	    _["row_idx"] = row_idx,
	    _["x"] = x,
	    _["icov"] = tmp_icov_p
	  );
}


double thresholdl1(double x, double thr) {
  if (x > thr)
    return x - thr;
  else if (x < -thr)
    return x + thr;
  else
    return 0;
}
