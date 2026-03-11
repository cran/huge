// Rcpp thin wrapper for MB graph estimation — delegates to huge::mb() / huge::mb_scr()
#include <Rcpp.h>
#include "huge/huge_core.h"
using namespace Rcpp;

// Helper: convert core ColResult vector to R sparse-like output
static List mb_result_to_list(const huge::MBResult& res, int d) {
    int total_nnz = 0;
    for (int m = 0; m < d; m++) total_nnz += res.columns[m].vals.size();

    NumericVector x(total_nnz);
    IntegerVector col_cnz(d + 1), row_idx(total_nnz);
    col_cnz[0] = 0;
    int cnz = 0;
    for (int m = 0; m < d; m++) {
        const huge::ColResult& col = res.columns[m];
        for (size_t j = 0; j < col.vals.size(); j++) {
            x[cnz] = col.vals[j];
            row_idx[cnz] = col.indices[j];
            cnz++;
        }
        col_cnz[m + 1] = cnz;
    }
    return List::create(
        _["col_cnz"] = col_cnz,
        _["row_idx"] = row_idx,
        _["x"] = x
    );
}

//[[Rcpp::export]]
List SPMBscr(NumericMatrix S, NumericVector lambda, int nlambda, int d, int maxdf, IntegerMatrix idx_scr, int nscr)
{
    if (d <= 0 || nlambda <= 0 || maxdf <= 0 || nscr < 0) {
        int d_safe = d > 0 ? d : 0;
        return List::create(
            _["col_cnz"] = IntegerVector(d_safe + 1),
            _["row_idx"] = IntegerVector(0),
            _["x"] = NumericVector(0)
        );
    }

    huge::MBResult res = huge::mb_scr(S.begin(), d, lambda.begin(), nlambda,
                                      idx_scr.begin(), nscr);
    return mb_result_to_list(res, d);
}

//[[Rcpp::export]]
List SPMBgraph(NumericMatrix S, NumericVector lambda, int nlambda, int d, int maxdf)
{
    if (d <= 0 || nlambda <= 0 || maxdf <= 0) {
        int d_safe = d > 0 ? d : 0;
        return List::create(
            _["col_cnz"] = IntegerVector(d_safe + 1),
            _["row_idx"] = IntegerVector(0),
            _["x"] = NumericVector(0)
        );
    }

    huge::MBResult res = huge::mb(S.begin(), d, lambda.begin(), nlambda);
    return mb_result_to_list(res, d);
}
