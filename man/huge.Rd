\name{huge}
\alias{huge}
\title{
High-dimensional undirected graph estimation in one-step mode
}
\description{
The main function for high-dimensional undirected graph estimation. It allows the user to load \code{huge.npn(), huge.scr(),huge.subgraph()} sequentially as a pipeline to analyze data.
}
\usage{
huge(L, ind.group = NULL, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.1, 
alpha = 1, sym = "or", npn = TRUE, npn.func = "shrinkage", npn.thresh = NULL, 
approx = FALSE, scr = TRUE, scr.num = NULL, verbose = TRUE)
}
\arguments{
  \item{L}{
There are two options for input \code{L}: (1) An \code{n} by \code{d} data matrix \code{L} representing \code{n} observations in \code{d} dimensions. (2) A list \code{L} containing \code{L$data} as an \code{n} by \code{d} data matrix. The list \code{L} can also contain \code{L$theta} as the true graph adjacency matrix, please refer to the returned values for more details.
}
  \item{ind.group}{
A length \code{k} vector indexing a subset of all \code{d} variables. ONLY applicable when estimating a subgraph of the whole graph. The default value is \code{c(1:d)}.
}
  \item{lambda}{
A sequence of decresing positive numbers to control the regularization in Meinshausen & Buhlmann Graph Estimation via Lasso (GEL) when \code{approx = FALSE} or Graph Approximation via Correlation Thresholding (GACT) when {approx = TRUE}. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. When \code{approx = FALSE}, use with care - it is better to supply a decreasing sequence values than a single (small) value.
}
  \item{nlambda}{
The number of regularization/thresholding paramters. The default value is \code{30} if \code{approx = TRUE} and \code{10} if \code{approx = FALSE}.
}
  \item{lambda.min.ratio}{
The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.1}. 
}
  \item{alpha}{
The tuning parameter for the elastic-net regression. The default value is \code{1} (lasso). When some dense pattern exists in the graph or some variables are highly correlated, the elastic-net is encouraged for its grouping effect. ONLY applicable when \code{approx = FALSE}.
}
  \item{sym}{
Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}. ONLY applicable when \code{approx = FALSE}.
}
  \item{npn}{
If \code{npn = TRUE}, the nonparanormal transformation is applied to the input data \code{L} or \code{L$data}. The default value is \code{TRUE}.
}
  \item{npn.func}{
The transformation function used in the NonparaNormal(NPN) transformation. If \code{npn.func = "truncation"}, the truncated ECDF is applied. If \code{npn.func = "shrinkage"}, the shrunken ECDF is applied. The default value is \code{"shrinkage"}.ONLY applicable when \code{npn = TRUE}.
}
  \item{npn.thresh}{
The truncation threshold used in NPN transformation, ONLY applicable when \code{npn.func = "truncation"}. The default value is \cr \code{1/(4*(n^0.25)*sqrt(pi*log(n)))}.
}
  \item{approx}{
If \code{approx = FALSE}, GEL is implemented. If \code{approx = TRUE}, GACT is implemented. The defaulty value is \code{approx = FALSE}. 
}
  \item{scr}{
If \code{scr = TRUE}, the Graph Sure Screening(GSS) is applied to preselect the neighborhood before GEL. The default value is \code{TRUE} for \code{n<d} and \code{FALSE} for \code{n>=d}. ONLY applicable when \code{approx = FALSE}.
}
  \item{scr.num}{
The neighborhood size after the GSS (the number of remaining neighbors per node). ONLY applicable when \code{scr = TRUE}. The default value is \code{n-1}. An alternative value is \code{n/log(n)}. ONLY applicable when \code{scr = TRUE} and \code{approx = FALSE}.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
This function provides a general framework for high-dimensional undirected graph estimation. The package integrates data preprocessing (Gaussianization), neighborhood preselection, graph estimation, and model selection techniques into a pipeline. The NPN transformation is applied to preprocess the data and helps relax the normality assumption. The GSS subroutine preselects the graph neighborhood of each variable. In the graph estimation stage, the structure of either the whole graph or a pre-specified sub-graph is estimated by the GEL on the pre-screened data. In the case \code{d >> n or d >>k}, the computation is memory optimized and is targeted on larger-sclae problems (with d>3000). We also provide another efficient method, the GACT.
}
\value{
An object with S3 class \code{"huge"} is returned:  
  \item{data}{
The \code{n} by \code{d} data matrix from the input
}
  \item{theta}{
The true graph structure from the input. ONLY applicable when the input list L contains L$theta as the true graph structure.
}
  \item{ind.group}{
The \code{ind.group} from the input
}
  \item{ind.mat}{
The \code{scr.num} by \code{k} matrix with each column correspondsing to a variable in \code{ind.group} and contains the indices of the remaining neighbors after the GSS. ONLY applicable when \code{scr = TRUE} and \code{approx = FALSE}
}
  \item{lambda}{
The sequence of regularization parameters used in GEL or thresholding parameters in GACT.
}
  \item{alpha}{
The \code{alpha} from the input. ONLY applicable when \code{approx = FALSE}.
}
  \item{sym}{
The \code{sym} from the input. ONLY applicable when \code{approx = FALSE}.
}
  \item{npn}{
The \code{npn} from the input.
}
  \item{scr}{
The \code{scr} from the input. ONLY applicable when \code{approx = FALSE}.
}
  \item{graph}{
return "subgraph path" when \code{k<d} and "fullgraph path" when \code{k==d}.
}
  \item{path}{
A list of \code{k} by \code{k} adjacency matrices of estimated graphs is returned as the solution path corresponding to \code{lambda}.
}
  \item{sparsity}{
The sparsity levels of the solution path.
}
  \item{approx}{
The correlation graph estimation indicator from the input
}
  \item{rss}{
A \code{k} by \code{nlambda} matrix. Each row is corresponding to a variable in \code{ind.group} and contains all RSS's (Residual Sum of Squares) along the lasso solution path. ONLY applicable when \code{approx = FALSE}.
}
  \item{df}{
A \code{k} by \code{nlambda} matrix. Each row corresponds to a variable in \code{ind.group} and contains the number of nonzero coefficients along the lasso solution path. ONLY applicable when \code{approx = FALSE}.
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}
\references{
Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010\cr
Han Liu, John Lafferty and Larry Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research} (JMLR), Vol.10, Page 2295-2328, 2009 \cr
Jianqing Fan and Jinchi Lv. Sure independence screening for ultra-high dimensional feature space (with discussion). \emph{Journal of Royal Statistical Society B}, Vol.70, Page 849-911, 2008.\cr
Jerome Friedman, Trevor Hastie and Rob Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, Vol.33, No.1, 2008. \cr
Nicaolai Meinshausen and Peter Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics},  Vol.34, Page 1436-1462, 2006.
}

\note{
This function ONLY estimates the solution path. For more information about the optimal graph selection, please refer to \code{\link{huge.select}}.\cr
This function can ONLY work under the setting \code{d > 2} and \code{scr.num > 1}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{huge.generator}}, \code{\link{huge.npn}}, \code{\link{huge.scr}}, \code{\link{huge.subgraph}}, \code{\link{huge.select}}, \code{\link{huge.plot}}, \code{\link{huge.roc}}, \code{\link{lasso.stars}} and \code{\link{huge-package}}.
}

\examples{
#generate data
L = huge.generator(n = 200, d = 80, graph = "hub")
ind.group = c(1:50)

#subgraph solution path estimation with input as a list
out1 = huge(L, ind.group = ind.group)
summary(out1)
plot(out1)
plot(out1, align = TRUE)
huge.plot(out1$path[[3]])
plot(out1$lambda,out1$sparsity)

#subgraph solution path estimation using the GACT
out2 = huge(L$data, ind.group = ind.group, approx = TRUE)
summary(out2)
plot(out2)

#fullgraph solution path estimation using elastic net
out3 = huge(L, alpha = 0.3)
summary(out3)
plot(out3)
}