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
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty and Larry Wasserman \cr
Maintainer: Tuo Zhao<tzhao5@jhu.edu>;Han Liu<hanliu@cs.jhu.edu>
}

\references{
1.Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
2.Han Liu, John Lafferty and Larry Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research}, 2009
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