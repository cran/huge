\name{huge.roc}
\alias{huge.roc}
\title{
Draw ROC Curve for a solution path
}
\description{
Draws ROC curve for a solution path according to the true graph structure
}
\usage{
huge.roc(est, theta, ind.group, verbose = TRUE)
}
\arguments{
  \item{est}{
There are two choices for the input est: (1) A solution path. (2)An object with S3 class "huge" containing \code{est$path} as a solution path, \code{est$theta} as the true graph structure and \code{est$ind.group} containing indices from a subset of all variables.
}
  \item{theta}{
The true graph structure. Only applicable when the input \code{est} is a solution path.
}
  \item{ind.group}{
A \code{k} dimensional vector containing \code{k} indices from a subset of all variables. Only applicable when estimating a subgraph of the whole graph and the input \code{est} is a solution path. The default value is \code{c(1:d)}.

}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
For a lasso regression,  the number of nonzero coefficients is at most \code{n-1}. If the dimension \code{d>>n}, even when  lambda is small, the estimated graph may still be sparse. In this case, the AUC may not be appropriate to evaluate the performance. 
}
\value{
An object with S3 class "roc" is returned:
 \item{F1}{
The F1 scores along the solution path.
}
  \item{tp}{
The true positive rates along the solution path
}
  \item{fp}{
The false positive rates along the solution paths
}
  \item{AUC}{
Area under the ROC curve
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}
\references{
Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. Technical Report, Carnegie Mellon University, 2010
}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{huge}} and \code{\link{huge-package}}
}
\examples{
#generate data
L = huge.generator(d = 200, graph = "cluster", prob = 0.3)
ind.group = c(1:40)
out1 = huge(L, ind.group = ind.group)

#draw ROC curve
Z1 = huge.roc(out1)

#Maximum F1 score
max(Z1$F1)
}