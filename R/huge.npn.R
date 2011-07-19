#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.npn(): nonparanormal(npn) transofmration                         #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tzhao@jhu.edu> and <hanliu@cs.jhu.edu>                       #
# Date: Jul 15th 2011                                                   #
# Version: 1.1.0                                                        #
#-----------------------------------------------------------------------#

## Main function
huge.npn = function(x, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE){
	gcinfo(FALSE)
	n = nrow(x)
  	d = ncol(x)
  	x.col = colnames(x)
  	x.row = rownames(x)
  	
  	# Shrinkaage transformation
	if(npn.func == "shrinkage"){
		if(verbose) cat("Conducting the nonparanormal (npn) transformation via shrunkun ECDF....")
		
		x = qnorm(apply(x,2,rank)/(n+1))
		x = x/sd(x[,1])
		
		if(verbose) cat("done.\n")
		rm(n,d,verbose)
   		gc()	
	}
	
	# Truncation transformation
	if(npn.func == "truncation"){
		if(verbose) cat("Conducting nonparanormal (npn) transformation via truncated ECDF....")
		if(is.null(npn.thresh)) npn.thresh = 1/(4*(n^0.25)*sqrt(pi*log(n)))
		
		x = qnorm(pmin(pmax(apply(x,2,rank)/n, npn.thresh), 1-npn.thresh))
    	x = x/sd(x[,1])
    	
    	if(verbose) cat("done.\n")
    	rm(n,d,npn.thresh,verbose)
   		gc()
	}
	colnames(x) = x.col
	rownames(x) = x.row
	return(x)
}