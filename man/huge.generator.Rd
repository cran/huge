\name{huge.generator}
\alias{huge.generator}
\title{
Data generating function
}
\description{
Implements the data generator from multivariate normal distributions with different graph structures, including \code{"random"}, \code{"hub"}, \code{"clique"} or \code{"band"}.
}
\usage{
huge.generator(n = 200, d = 50, graph = "random", v = NULL, u = NULL, 
g = NULL, prob = 0.03, vis = FALSE, verbose = TRUE)
}
\arguments{
  \item{n}{
The number of observations (sample size). The default value is 200.
}
  \item{d}{
The number of variables (dimension). The default value is 50.
}
  \item{graph}{
The graph structure with 4 options: \code{"random"}, \code{"hub"}, \code{"clique"} and \code{"band"}.
}
  \item{v}{
The off-diagonal elements of the precision matrix, controlling the magnitude of partial correlations with \code{u}. The default value is \code{0.1} for \code{"band"} and  \code{0.3} for \code{"random"}, \code{"hub"}, \code{"cliques"}.
}
  \item{u}{
A positive number (with a default value 0.5) being added to the diagonal elements of the precision matrix, to control the magnitude of partial correlations. The default value is \code{0} for \code{"clique"} and  \code{1} for \code{"random"}, \code{"hub"}, \code{"band"}.
}
 \item{g}{
For \code{"clique"} or \code{"hub"} graph, \code{g} is the number of hubs or cliques in the graph: when \code{"graph = hub"} default value is about \code{d/8} if \code{d >= 16} and \code{2} if \code{d<16}; when \code{"graph = clique"} default value is about \code{d/5} if \code{d >= 10} and \code{2} if \code{d<10}. For \code{"band"} graph, \code{g} is the bandwidth and the default value is \code{1}. Not applicable to random graph.
}
 \item{prob}{
The probability that a pair of nodes has an edge. The default value is 0.03. Only applicable when the graph pattern is \code{"random"}.
}
  \item{vis}{
Visualize the adjacency matrix of the true graph structure, the graph pattern, the covariance matrix and the precision matrix. The default value is \code{FALSE}
}
\item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
Given the adjacency matrix \code{theta}, the graph patterns are generated as below:\cr\cr
(I) \code{"random"}: Each pair of off-diagonal elements are randomly set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} with probability \code{"prob"}, and \code{0} other wise. It results in about \code{d(d-1)prob/2} edges in the graph.\cr\cr
(II)\code{"hub"}:The row/columns are evenly partitioned into \code{g} disjoint groups. Each group is associated with a "center" row \code{i} in that group. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} if \code{j} also belongs to the same group as \code{i} and \code{0} otherwise. It results in \code{d - g} edges in the graph.\cr\cr
(III)"clique":The row/columns are evenly partitioned into \code{g} disjoint groups. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} if both \code{i} and \code{j} belong to the same group, and \code{0} other wise. It results in about \code{g(d/g)(d/g-1)/2} edges in the graph.\cr\cr
(IV)\code{"band"}: The off-diagonal elements are set to be \code{theta[i,j]=1} if \code{1<=|i-j|<=g} and \code{0} other wise. It results in \code{(2d-1-g)g/2} edges in the graph.\cr\cr
The adjacency matrix \code{theta} has all diagonal elements equal to \code{0}. To obtain a positive definite precision matrix, the smallest eigenvalue of \code{theta*v} is computed. Suppose \code{e} be the smallest eigenvalue and we let the precision matrix equals \code{theta*v+(|e|+0.01+t)I}. The covariance matrix is then computed to generate multivariate normal data.
}
\value{
An object with S3 class "sim" is returned:
\item{data}{
The \code{n} by \code{d} matrix for the generated data
}
\item{sigma}{
The covariance matrix for the generated data
}
  \item{omega}{
The precision matrix for the generated data
}
  \item{theta}{
The adjacency matrix of true graph structure (in sparse matrix representation) for the generated data
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@gmail.com>; Han Liu <hanliu@cs.jhu.edu>
}
\references{
Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
Jerome Friedman, Trevor Hastie and Rob Tibshiran. Applications of the lasso and grouped lasso to the estimation of sparse graphical models, \emph{Technical Report}, Stanford University, 2010
}
%%\Warning{
%%
%%}
\seealso{
\code{\link{huge}} and \code{\link{huge-package}}
}
\examples{
# band graph with bandwidth 3
L = huge.generator(n = 500, graph = "band", g = 3)
summary(L)
plot(L)

# random sparse graph
L = huge.generator(prob = 0.05, vis = TRUE)

# random dense graph
L = huge.generator(prob = 0.3, vis = TRUE)

# hub graph with 4 hubs
L = huge.generator(graph = "hub", g = 4, vis = TRUE)
}