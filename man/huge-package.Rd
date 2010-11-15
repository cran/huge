\name{huge-package}
\alias{huge-package}
\docType{package}
\title{
High-dimensional Undirected Graph Estimation
}
\description{
A package for high-dimensional undirected graph estimation
}
\details{
\tabular{ll}{
Package: \tab huge\cr
Type: \tab Package\cr
Version: \tab 0.8.1\cr
Date: \tab 2010-11-15\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The package "huge" provides 8 main functions:\cr
(1) The data generator creates random samples from multivariate normal distributions with different graph structures.\cr
(2) The nonparanormal transformation helps relax the normality assumption.\cr
(3) The graph screening helps reduce the computation burden by preselecting the neighborhood of each node by thresholding sample correlation.\cr
(4) Graph Approximation via Correlation Thresholding(GACT).\cr 
(5) Meinshausen & Buhlmann Graph Estimation via Lasso(GEL).\cr
(6) Model select using the Stability Approach to Regularization Selection(StARS)
(7) Model select using the Permutation Information Criterion (PIC).
(8) Model selection using the Extended Bayesian Information Criterion(EBIC)based on the pseudo-likelihood.\cr
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}
\references{
Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
Han Liu, John Lafferty and Larry Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research} (JMLR), Vol.10, Page 2295-2328, 2009 \cr
Jianqing Fan and Jinchi Lv. Sure independence screening for ultra-high dimensional feature space (with discussion). \emph{Journal of Royal Statistical Society B}, Vol.70, Page 849-91sss1, 2008.\cr
Jerome Friedman, Trevor Hastie and Rob Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, Vol.33, No.1, 2008. \cr
Jerome Friedman, Trevor Hastie and Rob Tibshiran. Applications of the lasso and grouped lasso to the estimation of sparse graphical models, \emph{Technical Report}, Stanford University, 2010.\cr
Nicaolai Meinshausen and Peter Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics},  Vol.34, Page 1436-1462, 2006.\cr
Han Liu, Kathryn Roeder and Larry Wasserman. Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
Jiahua Chen and Zehua Chen. Extended Bayesian information criterion for model selection with large model space. \emph{Biometrika}. 95, 759-771.\cr
}
\seealso{
\code{\link{huge.generator}}, \code{\link{huge.npn}}, \code{\link{huge.scr}}, \code{\link{huge.subgraph}}, \code{\link{huge}}, \code{\link{huge.plot}}, \code{\link{huge.roc}} and \code{\link{lasso.stars}}
}