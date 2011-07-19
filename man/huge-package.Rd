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
Version: \tab 1.1.0\cr
Date: \tab 2011-07-15\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The package "huge" provides 8 main functions:\cr
(1) the data generator creates random samples from multivariate normal distributions with different graph structures. Please refer to \code{\link{huge.generator}}.\cr
(2) the nonparanormal (npn) transformation helps relax the normality assumption. Please refer to \code{\link{huge.npn}}.\cr
(3) The correlation thresholding graph estimation. Please refer to \code{\link{huge}}.\cr 
(4) The Meinshausen-Buhlmann graph estimation. Please refer to \code{\link{huge}}.\cr
(5) The graphical Lasso algorithm using lossless screening rule. Please refer and \code{\link{huge}}.\cr
\cr
**Both (4) and (5) can be furether accelerated by the lossy screening rule preselecting the neighborhood of each node via thresholding sample correlation.
\cr
(6) The model selection using the stability approach to regularization selection. Please refer to \code{\link{huge.select}}.\cr
(7) The model selection using the rotation information criterion. Please refer to \code{\link{huge.select}}.\cr
(8) The model selection using the extended Bayesian information criterion. Please refer to \code{\link{huge.select}}.\cr
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tzhao5@jhu.edu>; Han Liu <hanliu@cs.jhu.edu>
}
\references{
1.Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
2.Han Liu, Kathryn Roeder and Larry Wasserman. Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
3.Jerome Friedman, Trevor Hastie and Rob Tibshirani. Applications of the lasso and grouped lasso to the estimation of sparse graphical models, \emph{Technical Report}, Stanford University, 2010.\cr
4.Han Liu, John Lafferty and Larry Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research}, 2009 \cr
5.Jianqing Fan and Jinchi Lv. Sure independence screening for ultra-high dimensional feature space (with discussion). \emph{Journal of Royal Statistical Society B}, 2008.\cr
6.Onureena Banerjee, Laurent El Ghaoui, Alexandre d'Aspremont: Model Selection Through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data. \emph{Journal of Machine Learning Research}, 2008.\cr
7.R. Foygel and M. Drton. Extended bayesian information criteria for gaussian graphical models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
8.Jerome Friedman, Trevor Hastie and Robert Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, 2008. \cr
9.Jerome Friedman, Trevor Hastie and Robert Tibshirani. Sparse inverse covariance estimation with the lasso, \emph{Biostatistics}, 2007.\cr
10.Nicolai Meinshausen and Peter Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics}, 2006.\cr
11. D. Witten and J. Friedman. A fast screening rule for the graphical lasso. \emph{Journal of Computational and Graphical Statistics}, to appear, 2011.
}
\seealso{
\code{\link{huge.generator}}, \code{\link{huge.npn}}, \code{\link{huge}}, \code{\link{huge.plot}} and \code{\link{huge.roc}}
}