\name{huge.npn}
\alias{huge.npn}

\title{
Nonparanormal(npn) transformation
}

\description{
Implements the Gausianization to help relax the assumption of normality.
}

\usage{
huge.npn(x, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
}

\arguments{
  \item{x}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
}
  \item{npn.func}{
The transformation function used in the npn transformation. If \code{npn.func = "truncation"}, the truncated ECDF is applied. If \code{npn.func = "shrinkage"}, the shrunken ECDF is applied. The default is \code{"shrinkage"}. If \code{npn.func = "skeptic"}, the nonparanormal skeptic is applied.
}
  \item{npn.thresh}{
The truncation threshold used in nonparanormal transformation, ONLY applicable when \code{npn.func = "truncation"}. The default value is \code{1/(4*(n^0.25)*} \code{sqrt(pi*log(n)))}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{
The nonparanormal extends Gaussian graphical models to semiparametric Gaussian copula models.Motivated by sparse additive models, the nonparanormal method estimates the Gasussian copula by marginally transforming the variables using smooth functions.Computationally, the estimation of a nonparanormal transformation is very efficient and only requires one pass of the data matrix.
}

\value{
\item{data}{
A \code{d} by \code{d} nonparanormal correlation matrix if \code{npn.func = "skeptic"}, and A \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} transformed dimensions other wise.
}
}

\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tzhao5@jhu.edu>
}

\references{
1.  T. Zhao and H. Liu. The huge Package for High-dimensional Undirected Graph Estimation in R. \emph{Journal of Machine Learning Research}, 2012\cr
2.  H. Liu, F. Han, M. Yuan, J. Lafferty and L. Wasserman. High Dimensional Semiparametric Gaussian Copula Graphical Models. \emph{Annals of Statistics},2012 \cr
3.  D. Witten and J. Friedman. New insights and faster computations for the graphical lasso. \emph{Journal of Computational and Graphical Statistics}, to appear, 2011.
4.  Han Liu, Kathryn Roeder and Larry Wasserman. Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
5.  R. Foygel and M. Drton. Extended bayesian information criteria for gaussian graphical models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
6.  H. Liu, J. Lafferty and L. Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research}, 2009 \cr
7.  J. Fan and J. Lv. Sure independence screening for ultra-high dimensional feature space (with discussion). \emph{Journal of Royal Statistical Society B}, 2008.\cr
8.  O. Banerjee, L. E. Ghaoui, A. d'Aspremont: Model Selection Through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data. \emph{Journal of Machine Learning Research}, 2008.\cr
9.  J. Friedman, T. Hastie and R. Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, 2008. \cr
10. J. Friedman, T. Hastie and R. Tibshirani. Sparse inverse covariance estimation with the lasso, \emph{Biostatistics}, 2007.\cr
11. N. Meinshausen and P. Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics}, 2006.\cr
}

\seealso{
  \code{\link{huge}} and \code{\link{huge-package}}.
}

\examples{
# generate nonparanormal data
L = huge.generator(graph = "cluster", g = 5)
L$data = L$data^5

# transform the data using the shrunken ECDF
Q = huge.npn(L$data)

# transform the non-Gaussian data using the truncated ECDF
Q = huge.npn(L$data, npn.func = "truncation")

# transform the non-Gaussian data using the truncated ECDF
Q = huge.npn(L$data, npn.func = "skeptic")

}