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
Version: \tab 1.0\cr
Date: \tab 2011-02-28\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The package "huge" provides 9 main functions:\cr
(1) the data generator creates random samples from multivariate normal distributions with different graph structures. Please refer to \code{\link{huge.generator}}.\cr
(2) the NonparaNormal (NPN) transformation helps relax the normality assumption. Please refer to \code{\link{huge.NPN}}.\cr
(3) The Graph SURE Screening (GSS) helps reduce the computation burden by preselecting the neighborhood of each node via thresholding sample correlation. Please refer to \code{\link{huge.MBGEL}}.\cr
(4) The Graph Estimation via Correlation Thresholding (GECT).Please refer to \code{\link{huge.GECT}} and \code{\link{huge}}.\cr 
(5) The Meinshausen & Buhlmann Graph Estimation via Lasso(MBGEL). Please refer to \code{\link{huge.MBGEL}} and \code{\link{huge}}.\cr
(6) The slightly modified Graphical Lasso (GLASSO) using sparse matrix representation. Please refer to \code{\link{huge.glassoM}} and \code{\link{huge}}.\cr
(7) The model selection using the Stability Approach to Regularization Selection(StARS). Please refer to \code{\link{huge.select}} and \code{\link{lasso.stars}}.\cr
(8) The model selection using the Rotation Information Criterion (RIC). Please refer to \code{\link{huge.select}}.\cr
(9) The model selection using the Extended Bayesian Information Criterion(EBIC). Please refer to \code{\link{huge.select}}.\cr
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}
\references{
1.Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
2.Han Liu, Kathryn Roeder and Larry Wasserman. Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models. \emph{Advances in Neural Information Processing Systems}(NIPS), 2010.\cr
3.Jerome Friedman, Trevor Hastie and Rob Tibshirani. Applications of the lasso and grouped lasso to the estimation of sparse graphical models, \emph{Technical Report}, Stanford University, 2010.\cr
4.Han Liu, John Lafferty and Larry Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research} (JMLR), 2009 \cr
5.Jianqing Fan and Jinchi Lv. Sure independence screening for ultra-high dimensional feature space (with discussion). \emph{Journal of Royal Statistical Society B}, 2008.\cr
6.Onureena Banerjee, Laurent El Ghaoui, Alexandre d'Aspremont: Model Selection Through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data. \emph{Journal of Machine Learning Research} (JMLR), 2008.\cr
7.Jiahua Chen and Zehua Chen. Extended Bayesian information criterion for model selection with large model space. \emph{Biometrika}, 2008.\cr
8.Jerome Friedman, Trevor Hastie and Robert Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, 2008. \cr
9.Jerome Friedman, Trevor Hastie and Robert Tibshirani. Sparse inverse covariance estimation with the lasso, \emph{Biostatistics}, 2007.\cr
10.Nicolai Meinshausen and Peter Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics}, 2006.
}
\seealso{
\code{\link{huge.generator}}, \code{\link{huge.NPN}}, \code{\link{huge.GECT}}, \code{\link{huge.MBGEL}}, \code{\link{huge.glassoM}}, \code{\link{huge}}, \code{\link{huge.plot}}, \code{\link{huge.roc}} and \code{\link{lasso.stars}}
}