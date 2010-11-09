#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.plot(): graph visualization function                             #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@gmail.com>; <hanliu@cs.jhu.edu>                     #
# Date: Nov 9th, 2010                                                  #
# Version: 0.7                                                          #
#-----------------------------------------------------------------------#

huge.plot = function(G, epsflag = FALSE, graph.name = "default", cur.num = 1, location=NULL){
	par(mfrow = c(1,1))
	G = as.matrix(G)
	if(sum(G)==0) cat("The graph is a null graph.")
	if(sum(G)>0){
	if(missing(location))	location = getwd()
	setwd(location)
	g = graph.adjacency(G, mode="undirected", diag=FALSE)
	layout.grid = layout.fruchterman.reingold(g)
	
   if(epsflag == TRUE)	postscript(paste(paste(graph.name, cur.num, sep=""), "eps", sep="."), width = 8.0, height = 8.0) 
   # the main plotting function from igraph            
	plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA)
	rm(g,location)	
   gc(gcinfo(verbose = FALSE))
   if(epsflag == TRUE) dev.off()
   }
}
