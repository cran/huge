#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.npn(): The nonparanormal transofmration                          #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@gmail.com>; <hanliu@cs.jhu.edu>                     #
# Date: Nov 9th, 2010                                                   #
# Version: 0.7                                                          #
#-----------------------------------------------------------------------#

# Main function
huge.npn = function(x, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE){
	
	n = nrow(x)
  	d = ncol(x)
  	xt = list()
  	xt$npn.func = npn.func
  	
  	if(verbose) cat("Conducting the nonparanormal", npn.func, "transformation....")
  	
  	# Shrinkaage transformation
	if(npn.func == "shrinkage"){
		
		# Normal scores shrinkage
		for (i in 1:d) x[,i] = qnorm(rank(x[,i])/(n + 1))
		
		# Standardization
		xt$data = x/sd(x[,1])
		
		if(verbose) cat("done.\n")
		rm(x, n, d, verbose)
   		gc(gcinfo(verbose = FALSE)) 	
	}
	
	# Truncation transformation
	if(npn.func == "truncation"){
		
		if(is.null(npn.thresh)) npn.thresh = 1/(4*(n^0.25)*sqrt(pi*log(n)))
		
		# Normal scores truncation
		for (i in 1:d) {
    		tmp = rank(x[,i])/n
    		tmp[tmp > (1-npn.thresh)]= 1-npn.thresh
    		tmp[tmp < npn.thresh ]   = npn.thresh
    		x[,i] = qnorm(tmp)
    	}
    	
    	# Standardization
    	xt$data = x/sd(x[,1])
    	
    	if(verbose) cat("done.\n")
    	
    	rm(n,d,tmp,npn.func,npn.thresh,verbose)
   		gc(gcinfo(verbose = FALSE))
	}
	
	# Output class "npn"
	class(xt) = "npn"
	return(xt)
}

# Default print function for class "npn" 
print.npn = function(x, ...){
	cat("Nonparanormal transformation type:", x$npn.func,"\n")
}

# Default summary function for class "npn"
summary.npn = function(object, ...){
	cat("Nonparanormal transformation type:", object$npn.func,"\n")
}

# Default plot function for class "npn"
plot.npn = function(x, ...){
	cat("No plot information available.")
}