\name{huge.select}
\alias{huge.select}

\title{
Model selection for high-dimensional undirected graph estimation
}

\description{
Implements the regularization parameter selection for high dimensional undirected graph estimation. The optional approaches are Permutational Information Criterion (PIC), Stability Approach to Regularization Selection (StARS) and Extended Bayesian Information Criterion(EBIC).
}

\usage{
huge.select(est, criterion = NULL, r.num = 200, EBIC.gamma = 0.5, 
stars.thresh = 0.1, stars.subsample.ratio = NULL, stars.rep.num = 20, 
verbose = TRUE)
}

\arguments{
  \item{est}{
An object with S3 class \code{"huge"} (output from \code{\link{huge}})
}
  \item{criterion}{
Model selection criterion. For Meinshausen & Buhlmann Graph Estimation via Lasso (GEL), all 3 options are available, \code{"PIC"}, \code{"EBIC"} and \code{"stars"}. For Graph Approximation via Correlation Thresholding (GACT), \code{"stars"} is the Only option. For Graphical Lasso (GLASSO), \code{"EBIC"} and \code{"stars"} are available. The default value is \code{"PIC"} for the GEL, \code{"stars"} for the GACT, \code{EBIC} for the GLASSO.
}
  \item{r.num}{
The number of random permutations in PIC selection. The default value is 200. ONLY applicable when \code{criterion = "PIC"}
}
  \item{EBIC.gamma}{
The tuning parameter for the EBIC. The default value is 0.5. Only applicable when \code{est$method = "GEL"} or \code{est$method = "GLASSO"} and \code{criterion = "EBIC"}.
}
  \item{stars.thresh}{
The variability threshold in the StARS. The default value is \code{0.1}. An alternative value is \code{0.05}. Only applicable when \code{criterion = "stars"}.
}
  \item{stars.subsample.ratio}{
The subsampling ratio. The default value is \code{10*sqrt(n)/n} when \code{n>144} and \code{0.8} when \code{n<=144}, where \code{n} is the sample size. Only applicable when \code{criterion = "stars"}.
}
  \item{stars.rep.num}{
The number of subsampling for the StARS. The default value is \code{20}.Only applicable when \code{criterion = "stars"}
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{
The StARS is a natural way to select optimal regularization parameter for all three estimation methods. It selects the optimal graph by variability of subsamplings and tends to overselect edges in Gaussian graphical models. Besides selecting the regularization parameters, the StARS can also provide an additional estimated graph by merging the corresponding subsampled graphs using the frequency counts. The subsampling procedure in StARS may NOT be very efficient, we also proved the recent developed highly efficient PIC approach. Instead of tuning over a grid by cross-validation or subsampling, we directly estimate the optimal regularization paramter based on random permutations. However, the PIC usually has very good empirical performances but suffers from underselections sometimes. Therefore, we suggest if user are sensitive of false negative rates, they should either consider increasing \code{r.num} or applying the StARS to model selection. The extended BIC is another competive approach, but the \code{EBIC.gamma} can only be tuned by experience. The EBIC score for the GEL is based on pseudo-likelihood and the theoretical properties have NOT been justified yet.
}

\value{
An object with S3 class "select" is returned:
  \item{refit}{
The optimal graph selected from the solution path
}
  \item{merge}{
The solution path estimated by merging the subsampling paths. Only applicable when the input \code{criterion = "stars"}.
}
  \item{variability}{
The variability along the subsampling paths. Only applicable when the input \code{criterion = "stars"}.
}
  \item{EBIC.scores}{
Extended BIC scores for regularization parameter selection. Only applicable when \code{criterion = "EBIC"}.
}
  \item{opt.index}{
The index of the selected regularization parameter. NOT applicable when the input \code{criterion = "PIC"}
}
  \item{opt.lambda}{
The selected regularization/thresholding parameter.
}
  \item{opt.sparsity}{
The sparsity level of \code{"refit"}.
}
  \item{graph}{
return "subgraph" when \code{k<d} and "fullgraph" when \code{k==d}
}

and anything else inluded in the input \code{est}

}

\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}

\references{
1. Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. Technical Report, Carnegie Mellon University, 2010\cr
2. Han Liu, Kathryn Roeder and Larry Wasserman. Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models. \emph{Advances in Neural Information Processing Systems} (NIPS), 2010.\cr
Jiahua Chen and Zehua Chen. Extended Bayesian information criterion for model selection with large model space. \emph{Biometrika}, 2008.
}

\seealso{
\code{\link{huge}} and \code{\link{huge-package}}.
}

\examples{
#generate data
L = huge.generator(d = 200, graph="hub")
out.GEL = huge(L)
out.GACT = huge(L, method = "GACT")
out.GLASSO = huge(L, method = "GLASSO")

#model selection using PIC
out.select = huge.select(out.GEL)
summary(out.select)
plot(out.select)

#model selection using stars
out.select = huge.select(out.GACT, stars.rep.num = 5)
summary(out.select)
plot(out.select)

#model selection using EBIC
out.select = huge.select(out.GLASSO)
summary(out.select)
plot(out.select)
}