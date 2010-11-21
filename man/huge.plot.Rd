\name{huge.plot}
\alias{huge.plot}

\title{
Graph visualization
}

\description{
Implements the graph visualization using adjacency matrix. It can automatic organize 2D embedding layout.
}

\usage{
huge.plot(G, epsflag = FALSE, graph.name = "default", cur.num = 1, 
location)
}

\arguments{
  \item{G}{
The adjaceny matrix corresponding to the graph.
}
  \item{epsflag}{
If \code{epsflag = TRUE}, save the plot as an eps file in the target directory. The default value is \code{FALSE}.
}
  \item{graph.name}{
The name of the output eps files. The default value is "default".
}
  \item{cur.num}{
The number of plots saved as eps files. Only applicale when \code{epsflag = TRUE}. The default value is 1.
}
  \item{location}{
Target directory. The default value is the current working directory.
}
}

\details{
The user can change \code{cur.num} to plot several figures and select the best one. The implementation is based on the popular package "igraph".
}

\author{
Tuo Zhao, Han Liu, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
Maintainers: Tuo Zhao<tourzhao@andrew.cmu.edu>; Han Liu <hanliu@cs.jhu.edu>
}

\references{
1. Tuo Zhao and Han Liu. HUGE: A Package for High-dimensional Undirected Graph Estimation. \emph{Technical Report}, Carnegie Mellon University, 2010
}

\seealso{
\code{\link{huge}} and \code{\link{huge-package}}
}

\examples{
## visualize the hub graph
L = huge.generator(graph = "hub")
huge.plot(L$theta)

## visualize the band graph
L = huge.generator(graph = "band",g=5)
huge.plot(L$theta)

## visualize the cluster graph
L = huge.generator(graph = "cluster")
huge.plot(L$theta)

#show working directory
getwd()
#plot 5 graphs and save the plots as eps files in the working directory  
huge.plot(L$theta, epsflag = TRUE, cur.num = 5)
}