// Rcpp thin wrapper for TIGER estimation — delegates to huge::tiger()
#include <Rcpp.h>
#include "huge/huge_core.h"
using namespace Rcpp;

//[[Rcpp::export]]
List SPMBgraphsqrt(NumericMatrix data, NumericVector lambda, int nlambda, int d)
{
    int n = data.nrow();
    if (n <= 0 || d <= 0 || nlambda <= 0) {
        int d_safe = d > 0 ? d : 0;
        return List::create(
            _["col_cnz"] = IntegerVector(d_safe + 1),
            _["row_idx"] = IntegerVector(0),
            _["x"] = NumericVector(0),
            _["icov"] = List(0)
        );
    }

    huge::TigerResult res = huge::tiger(data.begin(), n, d, lambda.begin(), nlambda);

    // Reuse same CSC conversion pattern as SPMBgraph.cpp
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

    // Convert icov matrices
    List icov_list(nlambda);
    const size_t mat_bytes = sizeof(double) * static_cast<size_t>(d) * d;
    for (int i = 0; i < nlambda; i++) {
        NumericMatrix ic(d, d);
        std::memcpy(ic.begin(), res.icov[i].v.data(), mat_bytes);
        icov_list[i] = ic;
    }

    return List::create(
        _["col_cnz"] = col_cnz,
        _["row_idx"] = row_idx,
        _["x"] = x,
        _["icov"] = icov_list
    );
}
