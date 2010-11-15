\name{lasso.stars}
\alias{lasso.stars}
\title{
StARS Regularization Parameter Selection for Lasso
}
\description{
Implements the Stability Approach to Regularization Selection (StARS) for Lasso
}
\usage{
lasso.stars(x, y, rep.num = 20, lambda = NULL, nlambda = 100, 
lambda.min.ratio = 0.001, stars.thresh = 0.1, sample.ratio = NULL, 
alpha = 1, verbose = TRUE)
}
\arguments{
  \item{x}{
The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
}
  \item{y}{
The \code{n}-dimensional response vector
}
  \item{rep.num}{
The number of subsampling for StARS. The default value is \code{20}.
}
  \item{lambda}{
A sequence of decresing positive numbers to control regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. Use with care - it is better to supply a decreasing sequence values than a single (small) value.
}
  \item{nlambda}{
The number of regularization paramters. The default value is \code{100}.
}
  \item{lambda.min.ratio}{
The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.001}. 
}
  \item{stars.thresh}{
The threshold of the variability in StARS. The default value is \code{0.1}. The alternative value is \code{0.05}. Only applicable when \code{criterion = "stars"}
}
  \item{sample.ratio}{
The subsampling ratio. The default value is \code{10*sqrt(n)/n} when \code{n>144} and \code{0.8} when \code{n<=144}, where \code{n} is the sample size.
}
  \item{alpha}{
The tuning parameter for the elastic-net regression. The default value is \code{1} (lasso).
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
StARS selects the optimal regularization parameter based on the variability of the solution path. It chooses the least sparse graph among all solutions with the same variability. An alternative threshold \code{0.05} is chosen under the assumption that the model is correctly specified. In applications, the model  is usually an approximation of the true model, \code{0.1} is a safer choice.
}
\value{
An object with S3 class "stars" is returned:
  \item{path}{
The solution path of regression coefficients (in an \code{d} by \code{nlambda} matrix)
}
  \item{lambda}{
The regularization parameters used in Lasso
}
  \item{opt.index}{
The index of the optimal regularization parameter.
}
  \item{opt.beta}{
The optimal regression coefficients.
}
  \item{opt.lambda}{
The optimal regularization parameter.
}
  \item{Variability}{
The variability along the solution path.
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}
\references{
Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. Technical Report, Carnegie Mellon University, 2010\cr
Han Liu, Kathryn Roeder and Larry Wasserman. Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
Jerome Friedman, Trevor Hastie and Rob Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, Vol.33, No.1, 2008.
}

\note{
This function can only work under the setting when \code{d>1}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{huge.select}}, \code{\link{glmnet}}, \code{\link{huge}} and \code{\link{huge-package}}
}
\examples{
#generate data
x = matrix(rnorm(50*80),50,80)
beta = c(3,2,1.5,rep(0,77))
y = rnorm(50) + x\%*\%beta

#StARS for Lasso
z1 = lasso.stars(x,y)
summary(z1)
plot(z1)

#StARS for Lasso
z2 = lasso.stars(x,y, stars.thresh = 0.05)
summary(z2)
plot(z2)

#StARS for Lasso
z3 = lasso.stars(x,y,rep.num = 50)
summary(z3)
plot(z3)
}
