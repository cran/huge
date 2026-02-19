#include "math.h"
#include <algorithm>
#include <limits>
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double RIC(NumericMatrix &X, int d, int n, NumericVector &r, int t)
{
    if(d <= 1 || n <= 0 || t <= 0 || r.size() == 0)
        return 0.0;
    int i,j,k,m;
    int tmp_r;
    int t_eff = std::min(t, static_cast<int>(r.size()));

    double lambda_min,lambda_max,tmp;

    lambda_min = std::numeric_limits<double>::infinity();

    for(i=0;i<t_eff;i++)
    {
        tmp_r = static_cast<int>(r[i]);
        if(tmp_r < 0) tmp_r = 0;
        if(tmp_r > n) tmp_r = n;
        int split = n - tmp_r;
        lambda_max = 0;
        for(j=0;j<d;j++) {
            for(k=j+1;k<d;k++) {
                tmp = 0;
                for(m=0;m<split;m++)
                    tmp = tmp + X(m+tmp_r, j)*X(m,k);
                for(m=split;m<n;m++)
                    tmp = tmp + X(m-split, j)*X(m,k);
                tmp = fabs(tmp);
                if(tmp>lambda_max)
                    lambda_max = tmp;
            }
        }
        if(lambda_max<lambda_min)
            lambda_min = lambda_max;
    }
    return lambda_min;
}
