\name{huge.MBGEL}
\alias{huge.MBGEL}

\title{
Meinshausen & Buhlmann Graph Estimation via Lasso
}

\description{
Implements Meinshausen & Buhlmann Graph Estimation via Lasso (MBGEL). It estimates the neighborhood of each variable by fitting a collection of Lasso regression problems.
}

\usage{
huge.MBGEL(x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, scr = NULL, 
scr.num = NULL, idx.mat = NULL, sym = "or", verbose = TRUE)
}

\arguments{
  \item{x}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
}
  \item{idx.mat}{
A \code{scr.num} by \code{d} matrix. Each column contains the indices of the preslected neighborhood. Typical usage is to leave the input \code{idx.mat = NULL} and have the program compute its own \code{idx.mat} matrix based on \code{scr.num} by the Graph Sure Screening (GSS). Users can also specify a sequence to override this. Use with care, the indices are 0-based.

}
  \item{lambda}{
A sequence of decresing positive numbers to control regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. Use with care - it is better to supply a decreasing sequence values than a single (small) value.
}
  \item{nlambda}{
The number of regularization paramters. The default value is \code{10}.
}
  \item{lambda.min.ratio}{
The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.1}. 
}
  \item{scr}{
If \code{scr = TRUE}, the Graph Sure Screening(GSS) is applied to preselect the neighborhood for MBGEL. The default value is \code{TRUE} for \code{n<d} and \code{FALSE} for \code{n>=d}.
}
  \item{scr.num}{
The neighborhood size after the GSS (the number of remaining neighbors per node). ONLY applicable when \code{scr = TRUE}. The default value is \code{n-1}. An alternative value is \code{n/log(n)}.
}
  \item{sym}{
Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as neighbors for each other. The default value is \code{"or"}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{The MBGEL simplifies the precision matrix estimation to fitting a collection of Lasso regression problems by using each variable as response and the others as predictors. Unlike Graphical Lasso (GLASSO), it cannot numerically estimate the precision matrix but ONLY the underlying structure. It is computationally easier and can be more flexible in high-dimensional settings. It can be furether accelerated by the Graph Sure Screening (GSS). The GSS is applied to preselect the neighborhood under ultrahigh-dimensional setting. With the dimensionality reduced from ultra-high to a medium level (usually below the sample size), the GSS can greatly reduce the computational burden and often achieves equally or better estimation without using the GSS. The implementation is based on C.
}

\value{
An object with S3 class \code{"MBGEL"} is returned:
  \item{path}{
A list of \code{k} by \code{k} adjacency matrices (in sparse matrix representation) of estimated graphs as the solution path corresponding to \code{lambda}.
}
  \item{lambda}{
The sequence of regularization parameters used the graph estimation.
}
  \item{rss}{
A \code{k} by \code{nlambda} matrix. Each row is corresponding to a variable in \code{ind.group} and contains all RSS's (Residual Sum of Squares) along the lasso solution path.
}
  \item{df}{
A \code{k} by \code{nlambda} matrix. Each row corresponds to a variable in \code{ind.group} and contains the number of nonzero coefficients along the lasso solution path.
}
  \item{sparsity}{
The sparsity levels of the graph path.
}
}

\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}

\references{
1. Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
2.Jianqing Fan and Jinchi Lv. Sure independence screening for ultra-high dimensional feature space (with discussion). \emph{Journal of Royal Statistical Society B}, 2008.\cr
3. Jerome Friedman, Trevor Hastie and Rob Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, 2008. \cr
4. Nicaolai Meinshausen and Peter Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics}, 2006.
}

\note{
This function can ONLY work under the setting \code{d > 2} and \code{is.matrix(idx.mat) == TRUE}.
}

\seealso{
\code{\link{huge}} and \code{\link{huge-package}}
}

\examples{
#generate data
L = huge.generator(n = 100, d = 200, graph = "hub")

#graph path estimation with the GSS
out = huge.MBGEL(L$data)
summary(out)
plot(out)

#graph path estimation with specified lambda.min.ratio and nlambda
out = huge.MBGEL(L$data, nlambda = 8, lambda.min.ratio = 0.05)
summary(out)
plot(out)

#graph path estimation without the GSS
sub.path = huge.MBGEL(L$data, scr = FALSE)
summary(sub.path)
plot(sub.path)
}