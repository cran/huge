\name{plot.huge}
\alias{plot.huge}
\title{
Plot function for S3 class "huge"
}
\description{
Plot sparsity level information and 3 typical sparse graphs from the graph path
}
\usage{
\method{plot}{huge}(x, align = FALSE, ...)
}
\arguments{
  \item{x}{
An object with S3 class \code{"huge"} 
}
  \item{align}{
If \code{align = FALSE}, 3 plotted graphs are aligned
}
  \item{\dots}{
System reserved (No specific usage)
}
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
\code{\link{huge}}
}