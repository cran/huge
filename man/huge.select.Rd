\name{huge.select}
\alias{huge.select}
\title{
Model selection for high-dimensional undirected graph estimation
}
\description{
Implements the regularization parameter selection for high dimensional undirected graph estimation. The optional approaches are stability selection (StARS) and a variant of the extended BIC.
}
\usage{
huge.select(est, criterion = NULL, EBIC.gamma = 0.5, stars.thresh = 0.1, 
sample.ratio = NULL, rep.num = 20, verbose = TRUE)
}
\arguments{
  \item{est}{
An object with S3 class \code{"huge"} (output from \code{\link{huge}})
}
  \item{criterion}{
Model selection criterion. If \code{criterion = "EBIC"}, the extended BIC is applied (only applicable when \code{"est$approx = FALSE"}). If \code{criterion = "stars"}, StARS is applied. The default value is \code{"EBIC"} when \code{est$approx = FALSE} and \code{"stars"} when \code{est$approx = TRUE} .
}
  \item{EBIC.gamma}{
The tuning parameter for the extended BIC criterion. The default value is 0.5. Only applicable when \code{est$approx = FALSE} and \code{criterion = "EBIC"}.
}
  \item{stars.thresh}{
The variability threshold in StARS selection. The default value is \code{0.1}. An alternative value is \code{0.05}. Only applicable when \code{criterion = "stars"}.
}
  \item{sample.ratio}{
The subsampling ratio. The default value is \code{10*sqrt(n)/n} when \code{n>144} and \code{0.8} when \code{n<=144}, where \code{n} is the sample size. Only applicable when \code{criterion = "stars"}.
}
  \item{rep.num}{
The number of subsampling for StARS selection. The default value is \code{20}.Only applicable when \code{criterion = "stars"}
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
The StARS is a natural way to select optimal regularization parameter for high dimensional undirected graphical models. It also provides an additional estimated graph by merging the corresponding subsampled graphs using the frequency counts. The subsampling procedure in StARS may not be very efficient, therefore we also provide another method, the extended BIC score based on pseudo-likelihood. However its theoretical properties have not been justified yet.
}
\value{
An object with S3 class "select" is returned:
  \item{refit}{
The optimal graph selected from the solution path
}
  \item{merge}{
The merged by estimates from different subsamples. Only applicable when the input \code{criterion = "stars"}.
}
  \item{EBIC.scores}{
Extended BIC scores for regularization parameter selection. Only applicable when \code{criterion = "EBIC"}.
}
  \item{opt.index}{
The index of the selected regularization parameter.
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
Maintainers: Tuo Zhao<tourzhao@gmail.com>; Han Liu <hanliu@cs.jhu.edu>
}
\references{
Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. Technical Report, Carnegie Mellon University, 2010\cr
Han Liu, Kathryn Roeder and Larry Wasserman. Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
Jiahua Chen and Zehua Chen. Extended Bayesian information criterion for model selection with large model space. \emph{Biometrika}. 95, 759-771.
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{huge.stars}}, \code{\link{huge}} and \code{\link{huge-package}}.
}
\examples{
#generate data
L = huge.generator(graph="hub")

#subset indices
ind.group = c(1:30)

#estimate subgraph solution path using Meinshausen & Buhlmann graph estimation
out.huge = huge(L,ind.group = ind.group)


#model selection using extended BIC scores
out.select = huge.select(out.huge)
summary(out.select)
plot(out.select)

#model selection using stars
out.select = huge.select(out.huge, criterion = "stars", rep.num = 5)
summary(out.select)
plot(out.select)

#estimate subgraph solution path using GECA
out.approx = huge(L,ind.group = ind.group, approx = TRUE)

#model selection using stars
out.select = huge.select(out.approx, rep.num = 10)
summary(out.select)
plot(out.select)
}