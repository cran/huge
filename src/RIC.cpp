// Rcpp thin wrapper for RIC — delegates to huge::ric()
#include <Rcpp.h>
#include "huge/huge_core.h"
using namespace Rcpp;

//[[Rcpp::export]]
double RIC(NumericMatrix &X, int d, int n, NumericVector &r, int t)
{
    if (d <= 1 || n <= 0 || t <= 0 || r.size() == 0)
        return 0.0;

    int t_eff = std::min(t, static_cast<int>(r.size()));

    // Convert r from double to int (R passes as numeric)
    std::vector<int> r_int(t_eff);
    for (int i = 0; i < t_eff; i++)
        r_int[i] = static_cast<int>(r[i]);

    // X is column-major in R (NumericMatrix)
    return huge::ric(X.begin(), n, d, r_int.data(), t_eff);
}
