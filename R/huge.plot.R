#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.plot(): graph visualization                                      #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 12th, 2010                                                  #
# Version: 0.8                                                          #
#-----------------------------------------------------------------------#

huge.plot = function(G, epsflag = FALSE, graph.name = "default", cur.num = 1, location=NULL){
	if(missing(location))	location = getwd()
	setwd(location)
	g = graph.adjacency(as.matrix(G), mode="undirected", diag=FALSE)
	layout.grid = layout.fruchterman.reingold(g)
	
   	if(epsflag == TRUE)	postscript(paste(paste(graph.name, cur.num, sep=""), "eps", sep="."), width = 8.0, height = 8.0)             
	par(mfrow = c(1,1))
	plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA)
	rm(g,location)	
   	gc(gcinfo(verbose = FALSE))
   	if(epsflag == TRUE) dev.off()
}
