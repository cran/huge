\name{huge.npn}
\alias{huge.npn}
\title{
Nonparanormal transformation
}
\description{
Implements the nonparanormal transformation to relax the normality assumption.
}
\usage{
huge.npn(x, npn.func = "shrinkage", npn.thresh, verbose = TRUE)
}
\arguments{
  \item{x}{
The \code{n} by \code{d} data matrix representing n observations in d dimensions
}
  \item{npn.func}{
The transformation function used in the nonparanormal transformation. If \code{npn.func = "truncation"}, the truncated ECDF is applied. If \code{npn.func = "shrinkage"}, the shrunken ECDF is applied. The default is \code{"shrinkage"}.
}
  \item{npn.thresh}{
The truncation threshold used in nonparanormal transformation, only applicable when \code{npn.func = "truncation"}. The default value is \code{1/(4*(n^0.25)*} \code{sqrt(pi*log(n)))}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
The transformed data are already standardized as sample mean zero and unit variance.
}
\value{
An object with S3 class "npn" is returned:
  \item{data}{
The \code{n} by \code{d} data matrix representing n observations in d transformed dimensions
}
  \item{npn.func}{
The \code{npn.func} from the input
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty and Larry Wasserman \cr
Maintainer: Tuo Zhao<tourzhao@gmail.com>;Han Liu<hanliu@cs.jhu.edu>
}
\references{
Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
Han Liu, John Lafferty and Larry Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research} (JMLR), Vol.10, Page 2295-2328, 2009
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{huge}} and \code{\link{huge-package}}.
}

\examples{
n = 100
L = huge.generator(n=n, graph = "hub")

# nonparanormal transformation using the shrunken ECDF
Q = huge.npn(L$data)

# nonparanormal transformation using the truncated ECDF
Q = huge.npn(L$data, npn.func = "truncation")


# nonparanormal transformation using truncated ECDF with specific threhold
Q = huge.npn(L$data, npn.func = "truncation", npn.thresh = 1/n)
}