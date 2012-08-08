\name{huge.roc}
\alias{huge.roc}

\title{
Draw ROC Curve for a graph path
}

\description{
Draws ROC curve for a graph path according to the true graph structure
}

\usage{
huge.roc(path, theta, verbose = TRUE)
}

\arguments{
  \item{path}{
A graph path.
}
  \item{theta}{
The true graph structure.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{
To avoid the horizontal oscillation, false positive rates is automaically sorted in the ascent oder and true positive rates also follow the same order.
}

\note{
For a lasso regression,  the number of nonzero coefficients is at most \code{n-1}. If \code{d>>n}, even when regularization parameter is very small, the estimated graph may still be sparse. In this case, the AUC may not be a good choice to evaluate the performance.
}

\value{
An object with S3 class "roc" is returned:
 \item{F1}{
The F1 scores along the graph path.
}
  \item{tp}{
The true positive rates along the graph path
}
  \item{fp}{
The false positive rates along the graph paths
}
  \item{AUC}{
Area under the ROC curve
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
\code{\link{huge}} and \code{\link{huge-package}}
}

\examples{
#generate data
L = huge.generator(d = 200, graph = "cluster", prob = 0.3)
out1 = huge(L$data)

#draw ROC curve
Z1 = huge.roc(out1$path,L$theta)

#Maximum F1 score
max(Z1$F1)
}