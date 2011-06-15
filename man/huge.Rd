\name{huge}
\alias{huge}

\title{
High-dimensional undirected graph estimation
}

\description{
The main function for high-dimensional undirected graph estimation. Three graph estimation methods, including Meinshausen & Buhlmann Graph Estimation via Lasso (MBGEL), Graphical Lasso (GLASSO) and Graph Estimation via Correlation Thresholding (GECT), are available for data analysis.
}

\usage{
huge(x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, method = "mbgel", scr = NULL, scr.num = NULL, cov.output = FALSE, sym = "or", verbose = TRUE)
}

\arguments{
  \item{x}{
There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).

}
  \item{lambda}{
A sequence of decresing positive numbers to control the regularization in MBGEL and GLASSO, or the thresholding in GECT. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. When \code{method = "mbgel"} use with care - it is better to supply a decreasing sequence values than a single (small) value.
}
  \item{nlambda}{
The number of regularization/thresholding paramters. The default value is \code{30} if \code{method = "gect"} and \code{10} if \code{method = "mbgel"} or \code{method = "glasso"}.
}
  \item{lambda.min.ratio}{
If \code{method = "mbgel"} or \code{method = "glasso"}, it is the smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "gect"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mbgel"} or \code{method = "glasso"}, and 0.05 \code{method = "gect"}.
}
  \item{method}{
Graph estimation methods with 3 options: \code{"mbgel"}, \code{"gect"} and \code{"glasso"}. The defaulty value is \code{"mbgel"}. 
}
  \item{scr}{
If \code{scr = TRUE}, Graph Sure Screening (GSS) is applied to preselect the neighborhood before MBGEL. The default value is \code{TRUE} for \code{n<d} and \code{FALSE} for \code{n>=d}. ONLY applicable when \code{method = "mbgel"}.
}
  \item{scr.num}{
The neighborhood size after the GSS (the number of remaining neighbors per node). ONLY applicable when \code{scr = TRUE}. The default value is \code{n-1}. An alternative value is \code{n/log(n)}. ONLY applicable when \code{scr = TRUE} and \code{method = "mbgel"}.
}
  \item{cov.output}{
If \code{cov.output = TRUE}, the outpu will inlcude a path of estimated covariance matrices. ONLY applicable when \code{method = "glasso"}. Since the estimated covariance matrices are generally not sparse, please use it with care, or it may take much memory under high-dimensional setting. The default value is \code{FALSE}.
}
  \item{sym}{
Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}. ONLY applicable when \code{method = "mbgel"}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
The graph stucture is estimated using Meinshausen & Buhlmann Graph Estimation via Lasso (MBGEL) by default and it can be further accelerated by Graph SURE Screening (GSS) subroutine which preselects the graph neighborhood of each variable. In the case d >> n, the computation is memory optimized and is targeted on larger-scale problems. We also provide two alternative approaches for the graph estimation stage:(1) Graph Estimation via Correlation Thresholding (GECT) which is highly efficient and (2) A slightly modified Graphical Lasso (GLASSO) procedure in which the memory usage is optimized using sparse matrix output.
}
\value{
An object with S3 class \code{"huge"} is returned:  
  \item{data}{
The \code{n} by \code{d} data matrix or \code{d} by \code{d} sample covariance matrix from the input
}
  \item{cov.input}{
An indicator of the sample covariance. 
}
  \item{ind.mat}{
The \code{scr.num} by \code{k} matrix with each column correspondsing to a variable in \code{ind.group} and contains the indices of the remaining neighbors after the GSS. ONLY applicable when \code{scr = TRUE} and \code{approx = FALSE}
}
  \item{lambda}{
The sequence of regularization parameters used in MBGEL or thresholding parameters in GECT.
}
  \item{sym}{
The \code{sym} from the input. ONLY applicable when \code{method = "mbgel"}.
}
  \item{scr}{
The \code{scr} from the input. ONLY applicable when \code{method = "mbgel"}.
}
  \item{path}{
A list of \code{k} by \code{k} adjacency matrices of estimated graphs as a graph path corresponding to \code{lambda}.
}
  \item{sparsity}{
The sparsity levels of the graph path.
}
  \item{wi}{
A list of \code{d} by \code{d} precision matrices as an alternative graph path (numerical path) corresponding to \code{lambda}. ONLY applicable when {method = "glasso"}
}
  \item{w}{
A list of \code{d} by \code{d} estimated covariance matrices corresponding to \code{lambda}. ONLY applicable when \code{cov.output = TRUE} and {method = "glasso"}
}
  \item{method}{
The method used in the graph estimation stage.
}
  \item{df}{
If \code{method = "mbgel"}, it is a \code{k} by \code{nlambda} matrix. Each row contains the number of nonzero coefficients along the lasso solution path. If \code{method = "glasso"}, it is a \code{nlambda} dimensional vector containing the number of nonzero coefficients along the graph path \code{wi}.
}
  \item{loglik}{
A \code{nlambda} dimensional vector containing the likelihood scores along the graph path (\code{wi}). ONLY applicable when \code{method = "glasso"} 
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@gmail.com>; Han Liu <hanliu@cs.jhu.edu>
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
\code{\link{huge.generator}}, \code{\link{huge.select}}, \code{\link{huge.plot}}, \code{\link{huge.roc}}, and \code{\link{huge-package}}.
}

\examples{
#generate data
L = huge.generator(n = 200, d = 80, graph = "hub")

#graph path estimation using MBGEL
out1 = huge(L$data)
out1
plot(out1)				 #Not aligned	
plot(out1, align = TRUE) #Aligned
huge.plot(out1$path[[3]])

#graph path estimation using the sample covariance matrix as the input.
out1 = huge(cor(L$data))
out1
plot(out1)				 #Not aligned	
plot(out1, align = TRUE) #Aligned
huge.plot(out1$path[[3]])

#graph path estimation using GECT
out2 = huge(L$data,method = "gect")
out2
plot(out2)

#graph path estimation using GLASSO
out3 = huge(L$data, method = "glasso")
out3
plot(out3)
}