\name{huge.plot}
\alias{huge.plot}

\title{
Graph visualization
}

\description{
Implements the graph visualization using adjacency matrix. It can automatic organize 2D embedding layout.
}

\usage{
huge.plot(G, epsflag = FALSE, graph.name = "default", cur.num = 1, 
location)
}

\arguments{
  \item{G}{
The adjaceny matrix corresponding to the graph.
}
  \item{epsflag}{
If \code{epsflag = TRUE}, save the plot as an eps file in the target directory. The default value is \code{FALSE}.
}
  \item{graph.name}{
The name of the output eps files. The default value is "default".
}
  \item{cur.num}{
The number of plots saved as eps files. Only applicale when \code{epsflag = TRUE}. The default value is 1.
}
  \item{location}{
Target directory. The default value is the current working directory.
}
}

\details{
The user can change \code{cur.num} to plot several figures and select the best one. The implementation is based on the popular package "igraph".
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
\code{\link{huge}} and \code{\link{huge-package}}
}

\examples{
## visualize the hub graph
L = huge.generator(graph = "hub")
huge.plot(L$theta)

## visualize the band graph
L = huge.generator(graph = "band",g=5)
huge.plot(L$theta)

## visualize the cluster graph
L = huge.generator(graph = "cluster")
huge.plot(L$theta)

#show working directory
getwd()
#plot 5 graphs and save the plots as eps files in the working directory  
huge.plot(L$theta, epsflag = TRUE, cur.num = 5)
}