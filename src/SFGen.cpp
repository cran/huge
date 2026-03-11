// Rcpp thin wrapper for scale-free graph generator — delegates to huge::sfgen()
#include <Rcpp.h>
#include "huge/huge_core.h"
using namespace Rcpp;

//[[Rcpp::export]]
List SFGen(int d0, int d)
{
    if (d <= 0 || d0 <= 0 || d0 > d)
        stop("Invalid input: require 0 < d0 <= d.");

    IntegerMatrix G(d, d);

    // Generate random values using R's RNG
    int nrand = d - d0;
    std::vector<double> rands(nrand > 0 ? nrand : 0);

    GetRNGstate();
    for (int i = 0; i < nrand; i++)
        rands[i] = R::runif(0.0, 1.0);
    PutRNGstate();

    huge::sfgen(d0, d, G.begin(), rands.data());

    return List::create(_["G"] = G);
}
