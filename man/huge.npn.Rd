\name{huge.npn}
\alias{huge.npn}
\title{
NonparaNormal(NPN) transformation
}
\description{
Implements the NonparaNormal transformation (Gausianization) to relax the assumption of normality.
}
\usage{
huge.npn(x, npn.func = "shrinkage", npn.thresh, verbose = TRUE)
}
\arguments{
  \item{x}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
}
  \item{npn.func}{
The transformation function used in the NPN transformation. If \code{npn.func = "truncation"}, the truncated ECDF is applied. If \code{npn.func = "shrinkage"}, the shrunken ECDF is applied. The default is \code{"shrinkage"}.
}
  \item{npn.thresh}{
The truncation threshold used in NPN transformation, ONLY applicable when \code{npn.func = "truncation"}. The default value is \code{1/(4*(n^0.25)*} \code{sqrt(pi*log(n)))}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
The NPN is a very helpful tool to relax the normality assumption. Computationally, fitting a high dimensional NPN is no more difficult than estimating a multivariate Gaussian and one can also apply to other existing packages such as glasso. The transformed data are already standardized as sample mean zero and unit variance.
}
\value{
An object with S3 class "npn" is returned:
  \item{data}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} transformed dimensions
}
  \item{ntdata}{
The original data matrix before the NPN transformation
}
  \item{npn.func}{
The transformation function used in the NPN transformation
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty and Larry Wasserman \cr
Maintainer: Tuo Zhao<tourzhao@andrew.cmu.edu>;Han Liu<hanliu@cs.jhu.edu>
}
\references{
Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
Han Liu, John Lafferty and Larry Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research} (JMLR), Vol.10, Page 2295-2328, 2009
}
%%\note{
%%  ~~further notes~~
%%}

\seealso{
  \code{\link{huge}} and \code{\link{huge-package}}.
}

\examples{
L = huge.generator(graph = "cluster", g = 5)

# transform the non-Gaussian data using the shrunken ECDF
Q = huge.npn(L$data^5)
summary(Q)
plot(Q)

# transform the non-Gaussian data using the truncated ECDF
Q = huge.npn(5^(L$data), npn.func = "truncation")
plot(Q)
}