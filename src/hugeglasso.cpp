// Rcpp thin wrapper for graphical lasso — delegates to huge::glasso()
//
// Memory strategy: NumericMatrix is an SEXP wrapper whose .begin() points
// directly into R's heap.  By memcpy-ing core results straight into
// NumericMatrix we avoid the extra Eigen::MatrixXd intermediate that Rcpp
// would then have to serialize again — cutting one full d×d copy per matrix.
#include <Rcpp.h>
#include "huge/huge_core.h"
using namespace Rcpp;

//[[Rcpp::export]]
List hugeglasso(NumericMatrix S, NumericVector lambda, bool scr, bool verbose, bool cov_output)
{
    int d = S.nrow();
    int nlambda = lambda.size();

    if (verbose) {
        if (scr) Rcout << "Conducting the graphical lasso (glasso) with lossy screening...";
        else     Rcout << "Conducting the graphical lasso (glasso) with lossless screening...";
    }

    // S is column-major in R; .begin() gives a direct const-free pointer.
    huge::GlassoResult res = huge::glasso(S.begin(), d, lambda.begin(), nlambda, scr, cov_output);

    if (verbose) Rcout << "done.\n";

    // Scalar vectors — single allocation each, no copies
    NumericVector loglik(nlambda), sparsity(nlambda);
    IntegerVector df(nlambda);
    for (int i = 0; i < nlambda; i++) {
        loglik[i]   = res.loglik[i];
        sparsity[i] = res.sparsity[i];
        df[i]       = res.df[i];
    }

    // Matrix outputs — write directly into R SEXP memory (1 memcpy each)
    List path(nlambda), icov(nlambda), cov;
    if (cov_output) cov = List(nlambda);
    const size_t mat_bytes = sizeof(double) * static_cast<size_t>(d) * d;
    for (int i = 0; i < nlambda; i++) {
        NumericMatrix p(d, d);
        std::memcpy(p.begin(), res.path[i].v.data(), mat_bytes);
        path[i] = p;

        NumericMatrix ic(d, d);
        std::memcpy(ic.begin(), res.icov[i].v.data(), mat_bytes);
        icov[i] = ic;

        if (cov_output) {
            NumericMatrix cv(d, d);
            std::memcpy(cv.begin(), res.cov[i].v.data(), mat_bytes);
            cov[i] = cv;
        }
    }

    List result;
    result["loglik"]   = loglik;
    result["sparsity"] = sparsity;
    result["df"]       = df;
    result["path"]     = path;
    result["icov"]     = icov;
    if (cov_output) result["cov"] = cov;
    return result;
}
