\name{huge.subgraph}
\alias{huge.subgraph}
\title{
Subgraph estimation using Meinshausen & Buhlmann Graph Estimation via Lasso
}
\description{
Implements Meinshausen & Buhlmann Graph Estimation via Lasso (GEL). It estimates the neighborhood of each variable by fitting a collection of lasso regression problems.
}
\usage{
huge.subgraph(x, ind.group = NULL, ind.mat = NULL, alpha = 1, lambda = NULL, 
n.lambda = 10, lambda.min = 0.1, sym = "or", verbose = TRUE)
}
\arguments{
  \item{x}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
}
  \item{ind.group}{
A \code{k} dimensional vector indexing a subset of all \code{d} variables. Only applicable when estimating a subgraph of the whole graph. The default value is \code{c(1:d)}.
}
  \item{ind.mat}{
A \code{scr.num} by \code{k} matrix. Each column corresponds to a variable in \code{ind.group} and contains the indices of the preslected neighborhood.
}
  \item{alpha}{
The tuning parameter for the elastic-net regression. The default value is \code{1} (lasso). When some dense pattern exists in the graph or some variables are highly correlated, the elastic-net are encouraged for its grouping effect.
}
  \item{lambda}{
A sequence of decresing positive numbers to control regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{n.lambda} and \code{lambda.min}. Users can also specify a sequence to override this. Use with care - it is better to supply a decreasing sequence values than a single (small) value.
}
  \item{n.lambda}{
The number of regularization paramters. The default value is \code{10}.
}
  \item{lambda.min}{
The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{n.lambda} starting from \code{MAX} to \code{lambda.min*MAX} in log scale. The default value is \code{0.1}. 
}
  \item{sym}{
Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected only when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as neighbors for each other. The default value is \code{"or"}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
Meinshausen & Buhlmann Graph Estimation via Laso simplifies the precision matrix estimation to fitting a collection of lasso regression problems. Unlike the graphical lasso it does not numerically estimate the precision matrix. It is computationally easier and can be more scalable: In many situation such as gene regulatory network analysis, we are only interested in the structure of the graph or a subgraph of the whole graph. 
}
\value{
An object with S3 class \code{"subgraph"} is returned:
  \item{path}{
A list of \code{k} by \code{k} adjacency matrices of estimated graphs is returned as the solution path corresponding to \code{lambda}.
}
  \item{lambda}{
The sequence of regularization parameters used the graph estimation.
}
  \item{rss}{
A \code{k} by \code{n.lambda} matrix. Each row is corresponding to a variable in \code{ind.group} and contains all RSS's (Residual Sum of Squares) along the lasso solution path.
}
  \item{df}{
A \code{k} by \code{n.lambda} matrix. Each row corresponds to a variable in \code{ind.group} and contains the number of nonzero coefficients along the lasso solution path.
}
  \item{sparsity}{
The sparsity levels of the solution path.
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@gmail.com>; Han Liu <hanliu@cs.jhu.edu>
}
\references{
Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. Technical Report, Carnegie Mellon University, 2010\cr
Jerome Friedman, Trevor Hastie and Rob Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, Vol.33, No.1, 2008. \cr
Nicaolai Meinshausen and Peter Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics},  Vol.34, Page 1436-1462, 2006.
}
\note{
This function can only work under the setting \code{d > 2} and \code{nrow(ind.mat) > 1}.
}
\seealso{
\code{\link{huge}} and \code{\link{huge-package}}
}
\examples{
#generate data
L = huge.generator(graph = "hub")

#subset indices
ind.group = c(1:40)

#generate the preselected neighborhood
ind.mat = huge.scr(L$data, ind.group)$ind.mat

#generated lambda
lambda = exp(seq(log(0.8),log(0.01),length=30))

#subgraph solution path estimation with the preselected neighborhood
G.path = huge.subgraph(L$data, ind.group = ind.group, ind.mat = ind.mat, lambda = lambda)

#subgraph solution path estimation with specified lambda.min and n.lambda
G.path = huge.subgraph(L$data, ind.group = ind.group, ind.mat = ind.mat, 
n.lambda = 20, lambda.min = 0.01)

#graph solution path estimation using elastic net
G.path = huge.subgraph(L$data, ind.group = ind.group, alpha = 0.7)
}