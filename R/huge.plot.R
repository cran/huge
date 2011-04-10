#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.plot(): graph visualization                                      #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Apr 10th 2011                                                   #
# Version: 1.0.1                                                        #
#-----------------------------------------------------------------------#

huge.plot = function(G, epsflag = FALSE, graph.name = "default", cur.num = 1, location=NULL){
	gcinfo(FALSE)
	if(missing(location))	location = getwd()
	setwd(location)
	g = graph.adjacency(as.matrix(G!=0), mode="undirected", diag=FALSE)
	layout.grid = layout.fruchterman.reingold(g)
	
   	if(epsflag == TRUE)	postscript(paste(paste(graph.name, cur.num, sep=""), "eps", sep="."), width = 8.0, height = 8.0)             
	par(mfrow = c(1,1))
	plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=2, vertex.label=NA)
	rm(g,location)	
   	gc()
   	if(epsflag == TRUE) dev.off()
}
