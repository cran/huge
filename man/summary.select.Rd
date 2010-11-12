\name{summary.select}
\alias{summary.select}
\title{
Summary function for S3 class "select"
}
\description{
Summary the information about the model usage, subgraph dimension, fullgraph dimension, model selection criterion, sparsity level of the optimal graph
}
\usage{
\method{summary}{select}(object, ...)
}
\arguments{
  \item{object}{
An object with S3 class \code{"select"}
}
  \item{\dots}{
System reserved (No specific usage)
}
}
\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}
\seealso{
\code{\link{huge.select}} and \code{\link{huge}}
}