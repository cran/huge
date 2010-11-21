\name{huge.subgraph}
\alias{huge.subgraph}

\title{
Subgraph estimation using Meinshausen & Buhlmann Graph Estimation via Lasso
}

\description{
Implements Meinshausen & Buhlmann Graph Estimation via Lasso (GEL). It estimates the neighborhood of each variable by fitting a collection of Lasso regression problems.
}

\usage{
huge.subgraph(x, ind.group = NULL, ind.mat = NULL, alpha = 1, lambda = NULL, 
nlambda = 10, lambda.min.ratio = 0.1, sym = "or", verbose = TRUE)
}

\arguments{
  \item{x}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
}
  \item{ind.group}{
A \code{k} dimensional vector indexing a subset of all \code{d} variables. ONLY applicable when estimating a subgraph of the whole graph. The default value is \code{c(1:d)}.
}
  \item{ind.mat}{
A \code{scr.num} by \code{k} matrix. Each column corresponds to a variable in \code{ind.group} and contains the indices of the preslected neighborhood.
}
  \item{alpha}{
The tuning parameter for the elastic-net regression. The default value is \code{1} (lasso). When some dense pattern exists in the graph or some variables are highly correlated, the elastic-net are encouraged for its grouping effect.
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
  \item{sym}{
Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as neighbors for each other. The default value is \code{"or"}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{
The GEL simplifies the precision matrix estimation to fitting a collection of Lasso regression problems by using each variable as response and the others as predictors. Unlike Graphical Lasso (GLASSO), it cannot numerically estimate the precision matrix but ONLY the underlying structure. It is computationally easier and can be more flexible in high-dimensional settings. In some situations such as gene regulatory network analysis, we are only interested in the structure of the graph or a subgraph of the whole graph and the GEL is more scalable than other existing algorithms. The implementation is based on the well-known highly computationally efficient package "glmnet".
}

\value{
An object with S3 class \code{"subgraph"} is returned:
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
The sparsity levels of the solution path.
}
}

\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}

\references{
1. Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
2. Jerome Friedman, Trevor Hastie and Rob Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, 2008. \cr
3. Nicaolai Meinshausen and Peter Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics}, 2006.
}

\note{
This function can ONLY work under the setting \code{d > 2} and \code{is.matrix(ind.mat) == TRUE}.
}

\seealso{
\code{\link{huge}} and \code{\link{huge-package}}
}

\examples{
#generate data
L = huge.generator(graph = "hub")
ind.group = c(1:30)
ind.mat = huge.scr(L$data, ind.group)$ind.mat
lambda = exp(seq(log(0.8),log(0.1),length=12))

#subgraph solution path estimation with the preselected neighborhood and specified lambda sequence
sub.path = huge.subgraph(L$data, ind.group = ind.group, ind.mat = ind.mat, lambda = lambda)
summary(sub.path)
plot(sub.path)

#subgraph solution path estimation with specified lambda.min.ratio and nlambda
sub.path = huge.subgraph(L$data, ind.group = ind.group, ind.mat = ind.mat, 
nlambda = 8, lambda.min.ratio = 0.01)
summary(sub.path)
plot(sub.path)

#graph solution path estimation using elastic net
sub.path = huge.subgraph(L$data, ind.group = ind.group, alpha = 0.7)
summary(sub.path)
plot(sub.path)
}