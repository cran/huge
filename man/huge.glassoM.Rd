\name{huge.glassoM}
\alias{huge.glassoM}

\title{
High-dimensional undirected graph estimation via Graphical Lasso
}

\description{
Estimates a graph path of sparse precision matrices using a sequence of regularization parameters.
}

\usage{
huge.glassoM(x, lambda = NULL, lambda.min.ratio = NULL, nlambda = NULL, cov.output = FALSE, verbose = TRUE)
}

\arguments{
  \item{x}{
There are 2 options:(1) The \code{n} by \code{d} data matrix. (2) A \code{d} by \code{d} sample covariance covariance matrix. The program automatically identifies the input matrix by checking the symmetry.(\code{n} is the sample size and \code{d} is the dimension).
}
  \item{lambda}{
A sequence of decresing positive numbers to control the regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this.}
  \item{nlambda}{
The number of regularization paramters. The default value is \code{10}.
}
  \item{lambda.min.ratio}{
The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.1}.
}
  \item{cov.output}{
If \code{cov.output = TRUE}, the outpu will inlcude a path of estimated covariance matrices. Since the estimated covariance matrices are generally not sparse, please use it with care, or it may take much memory under high-dimensional setting. The default value is \code{FALSE}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
The graphical lasso is a well-known computational efficient algorithm to estimate the sparse precision matrix:
\deqn{log(det(Omega)) + Tr(S*Omega) + \lambda |Omega|_1}
The implementation of this function is based on the popular package "glasso". Our \code{huge.glassoM()} is actually a modified version of Friedman's glassopath() in a memory efficient manner. Since in most case we are usually only interested in sparse results, instead of using a regular array, we store the graph path in a list with each element corresponding to an estimated precision matrix or graph in sparse matrix representation. This slightly change can help to reduce the memory usage in high dimensional cases. We also follow Dr. Friedman's warm start procedure and make this function computational efficient and stable.
}
\value{
An object with S3 class \code{"glassoM"} is returned:  
  \item{path}{
A list of \code{k} by \code{k} adjacency matrices of estimated graphs (in sparse matrix representation) as a solution path (graph path) corresponding to \code{lambda}.
}
  \item{cov.input}{
An indicator of the sample covariance. 
}
  \item{wi}{
A list of \code{d} by \code{d} estimated precision matrices (in sparse matrix representation) as an alternative graph path (numerical path) corresponding to \code{lambda}.
}
  \item{w}{
A list of \code{d} by \code{d} estimated covariance matrices corresponding to \code{lambda}. ONLY applicable when \code{cov.output = TRUE}
}
  \item{lambda}{
The sequence of regularization parameters used in GLASSO.
}
  \item{sparsity}{
The sparsity levels of the graph path (graph path).
}
  \item{loglik}{
A \code{nlambda} dimensional vector containing the likelihood scores along the graph path (\code{wi}).
}
  \item{df}{
A \code{nlambda} dimensional vector containing the number of nonzero coefficients along the graph path (\code{wi}).
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@gmail.com>; Han Liu <hanliu@cs.jhu.edu>
}

\references{
1.Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
2.Jerome Friedman, Trevor Hastie and Robert Tibshirani. Sparse inverse covariance estimation with the lasso. \emph{Biostatistics}, 2007.\cr
6.Onureena Banerjee, Laurent El Ghaoui, Alexandre d'Aspremont: Model Selection Through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data. \emph{Journal of Machine Learning Research} (JMLR), 2008.\cr
}

\seealso{
\code{\link{huge}}, \code{\link{glasso}}, \code{\link{glassopath}} and \code{\link{huge-package}}.
}

\examples{
#generate data
L = huge.generator(n = 200, d = 100, graph = "hub")

#graph graph path estimation
out1 = huge.glassoM(L$data)
plot(out1)

#Graph graph path estimation with specified nlambda and lambda.min.ratio
out3 = huge(L, method = "GLASSO", nlambda = 8, lambda.min.ratio = 0.2)
plot(out3)
}