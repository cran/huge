\name{huge.NPN}
\alias{huge.NPN}

\title{
NonparaNormal(NPN) transformation
}

\description{
Implements the NonparaNormal transformation (Gausianization) to help relax the assumption of normality.
}

\usage{
huge.NPN(x, NPN.func = "shrinkage", NPN.thresh, verbose = TRUE)
}

\arguments{
  \item{x}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
}
  \item{NPN.func}{
The transformation function used in the NPN transformation. If \code{NPN.func = "truncation"}, the truncated ECDF is applied. If \code{NPN.func = "shrinkage"}, the shrunken ECDF is applied. The default is \code{"shrinkage"}.
}
  \item{NPN.thresh}{
The truncation threshold used in NPN transformation, ONLY applicable when \code{NPN.func = "truncation"}. The default value is \code{1/(4*(n^0.25)*} \code{sqrt(pi*log(n)))}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{
The NPN is a very helpful tool to relax the normality assumption. It assumes the data come from a Gaussian cornopula and computationally, fitting a high dimensional NPN is no more difficult than estimating a multivariate Gaussian and one can also apply to other existing estimation algorithms. The output data are already standardized as sample mean zero and unit variance.
}

\value{
An object with S3 class "NPN" is returned:
  \item{data}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} transformed dimensions
}
  \item{ntdata}{
The original data matrix before the NPN transformation
}
  \item{NPN.func}{
The transformation function used in the NPN transformation
}
}

\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty and Larry Wasserman \cr
Maintainer: Tuo Zhao<tourzhao@andrew.cmu.edu>;Han Liu<hanliu@cs.jhu.edu>
}

\references{
1.Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
2.Han Liu, John Lafferty and Larry Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research} (JMLR), 2009
}

\seealso{
  \code{\link{huge}} and \code{\link{huge-package}}.
}

\examples{
# generate data
L = huge.generator(graph = "cluster", g = 5)

# transform the non-Gaussian data using the shrunken ECDF
Q = huge.NPN(L$data^5)
summary(Q)
plot(Q)

# transform the non-Gaussian data using the truncated ECDF
Q = huge.NPN(5^(L$data), NPN.func = "truncation")
plot(Q)
}