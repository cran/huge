\name{huge}
\alias{huge}

\title{
High-dimensional undirected graph estimation in one-step mode
}

\description{
The main function for high-dimensional undirected graph estimation. It allows the user to load \code{huge.NPN(), huge.GECT(), huge.MBGEL(), huge.glassoM()} sequentially as a pipeline to analyze data.
}

\usage{
huge(L, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, NPN = FALSE, NPN.func = "shrinkage", NPN.thresh = NULL, method = "MBGEL", 
scr = NULL, scr.num = NULL, cov.glasso = FALSE, sym = "or", verbose = TRUE)
}

\arguments{
  \item{L}{
There are two options for input \code{L}: (1) An \code{n} by \code{d} data matrix \code{L} representing \code{n} observations in \code{d} dimensions. (2) A list \code{L} containing \code{L$data} as an \code{n} by \code{d} data matrix. The list \code{L} can also contain \code{L$theta} as the true graph adjacency matrix, please refer to the returned values for more details.
}
  \item{lambda}{
A sequence of decresing positive numbers to control the regularization in Meinshausen & Buhlmann Graph Estimation via Lasso (MBGEL) and Graphical Lasso (GLASSO), or the thresholding in Graph Estimation via Correlation Thresholding (GECT). Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. When \code{method = "MBGEL"} use with care - it is better to supply a decreasing sequence values than a single (small) value.
}
  \item{nlambda}{
The number of regularization/thresholding paramters. The default value is \code{30} if \code{method = "GECT"} and \code{10} if \code{method = "MBGEL"} or \code{method = "GLASSO"}.
}
  \item{lambda.min.ratio}{
If \code{method = "MBGEL"} or \code{method = "GLASSO"}, it is the smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "GECT"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "MBGEL"} or \code{method = "GLASSO"}, and 0.05 \code{method = "GECT"}.
}
  \item{NPN}{
If \code{NPN = TRUE}, the nonparanormal transformation is applied to the input data \code{L} or \code{L$data}. The default value is \code{FALSE}.
}
  \item{NPN.func}{
The transformation function used in the NonparaNormal(NPN) transformation. If \code{NPN.func = "truncation"}, the truncated ECDF is applied. If \code{NPN.func = "shrinkage"}, the shrunken ECDF is applied. The default value is \code{"shrinkage"}. ONLY applicable when \code{NPN = TRUE}.
}
  \item{NPN.thresh}{
The truncation threshold used in NPN transformation, ONLY applicable when \code{NPN.func = "truncation"}. The default value is \cr \code{1/(4*(n^0.25)*sqrt(pi*log(n)))}.
}
  \item{method}{
Graph estimation methods with 3 options: \code{"MBGEL"}, \code{"GECT"} and \code{"GLASSO"}. The defaulty value is \code{"MBGEL"}. 
}
  \item{scr}{
If \code{scr = TRUE}, the Graph Sure Screening(GSS) is applied to preselect the neighborhood before MBGEL. The default value is \code{TRUE} for \code{n<d} and \code{FALSE} for \code{n>=d}. ONLY applicable when \code{method = "MBGEL"}.
}
  \item{scr.num}{
The neighborhood size after the GSS (the number of remaining neighbors per node). ONLY applicable when \code{scr = TRUE}. The default value is \code{n-1}. An alternative value is \code{n/log(n)}. ONLY applicable when \code{scr = TRUE} and \code{method = "MBGEL"}.
}
  \item{cov.glasso}{
If \code{cov.glasso = TRUE}, the outpu will inlcude a path of estimated covariance matrices. ONLY applicable when \code{method = "GLASSO"}. Since the estimated covariance matrices are generally not sparse, please use it with care, or it may take much memory under high-dimensional setting. The default value is \code{FALSE}.
}
  \item{sym}{
Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}. ONLY applicable when \code{method = "MBGEL"}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
It provides a general framework for high-dimensional undirected graph estimation. It integrates data preprocessing (Gaussianization), neighborhood screening, graph estimation, and model selection techniques into a pipeline. In preprocessing stage, the NonparaNormal(NPN) transformation is applied to help relax the normality assumption. In the graph estimation stage, the graph structure is estimated by the Meinshausen & Buhlmann Graph Estimation via Lasso (MBGEL) by default and it can be further accelerated by the Graph SURE Screening (GSS) subroutine which preselects the graph neighborhood of each variable. In the case d >> n, the computation is memory optimized and is targeted on larger-scale problems (with d>10000). We also provide two alternative approaches for the graph estimation stage:(1) Graph Estimation via Correlation Thresholding (GECT) which is highly efficient and (2) A slightly modified Graphical Lasso (GLASSO) procedure in which the memory usage is optimized using sparse matrix output.
}
\value{
An object with S3 class \code{"huge"} is returned:  
  \item{data}{
The \code{n} by \code{d} data matrix from the input
}
  \item{theta}{
The true graph structure from the input. ONLY applicable when the input list L contains L$theta as the true graph structure.
}
  \item{ind.mat}{
The \code{scr.num} by \code{k} matrix with each column correspondsing to a variable in \code{ind.group} and contains the indices of the remaining neighbors after the GSS. ONLY applicable when \code{scr = TRUE} and \code{approx = FALSE}
}
  \item{lambda}{
The sequence of regularization parameters used in MBGEL or thresholding parameters in GECT.
}
  \item{sym}{
The \code{sym} from the input. ONLY applicable when \code{method = "MBGEL"}.
}
  \item{NPN}{
The \code{NPN} from the input.
}
  \item{scr}{
The \code{scr} from the input. ONLY applicable when \code{method = "MBGEL"}.
}
  \item{path}{
A list of \code{k} by \code{k} adjacency matrices of estimated graphs as a graph path corresponding to \code{lambda}.
}
  \item{sparsity}{
The sparsity levels of the graph path.
}
  \item{wi}{
A list of \code{d} by \code{d} precision matrices as an alternative graph path (numerical path) corresponding to \code{lambda}. ONLY applicable when {method = "GLASSO"}
}
  \item{w}{
A list of \code{d} by \code{d} estimated covariance matrices corresponding to \code{lambda}. ONLY applicable when \code{cov.glasso = TRUE} and {method = "GLASSO"}
}
  \item{method}{
The method used in the graph estimation stage.
}
  \item{rss}{
A \code{k} by \code{nlambda} matrix. Each row is corresponding to a variable in \code{ind.group} and contains all RSS's (Residual Sum of Squares) along the lasso solution path. ONLY applicable when \code{method = "MBGEL"}.
}
  \item{df}{
If \code{method = "MBGEL"}, it is a \code{k} by \code{nlambda} matrix. Each row corresponds to a variable in \code{ind.group} and contains the number of nonzero coefficients along the lasso solution path. If \code{method = "GLASSO"}, it is a \code{nlambda} dimensional vector containing the number of nonzero coefficients along the graph path \code{wi}.
}
  \item{loglik}{
A \code{nlambda} dimensional vector containing the likelihood scores along the graph path (\code{wi}). ONLY applicable when
\code{method = "GLASSO"} 
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}

\references{
1.Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
2.Han Liu, John Lafferty and Larry Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research} (JMLR), 2009 \cr
3.Jianqing Fan and Jinchi Lv. Sure independence screening for ultra-high dimensional feature space (with discussion). \emph{Journal of Royal Statistical Society B}, 2008.\cr
4.Jerome Friedman, Trevor Hastie and Rob Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, 2008. \cr
5.Onureena Banerjee, Laurent El Ghaoui, Alexandre d'Aspremont: Model Selection Through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data. \emph{Journal of Machine Learning Research} (JMLR), 2008.\cr
6.Jerome Friedman, Trevor Hastie and Robert Tibshirani. Sparse inverse covariance estimation with the lasso. \emph{Biostatistics}, 2007.\cr
7.Nicolai Meinshausen and Peter Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics}, 2006.
}

\note{
This function ONLY estimates the graph path. For more information about the optimal graph selection, please refer to \code{\link{huge.select}}.\cr
}

\seealso{
\code{\link{huge.generator}}, \code{\link{huge.NPN}}, \code{\link{huge.GECT}}, \code{\link{huge.MBGEL}}, \code{\link{huge.glassoM}}, \code{\link{huge.select}}, \code{\link{huge.plot}}, \code{\link{huge.roc}}, \code{\link{lasso.stars}} and \code{\link{huge-package}}.
}

\examples{
#generate data
L = huge.generator(n = 200, d = 80, graph = "hub")

#graph path estimation with input as a list
out1 = huge(L)
summary(out1)
plot(out1)
plot(out1, align = TRUE)
huge.plot(out1$path[[3]])
plot(out1$lambda,out1$sparsity)

#graph path estimation using the GECT
out2 = huge(L$data,method = "GECT")
summary(out2)
plot(out2)

#graph path estimation using the GLASSO
out3 = huge(L, method = "GLASSO")
summary(out3)
plot(out3)
}