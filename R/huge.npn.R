#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.NPN(): NonparaNormal transofmration                              #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Feb 28th 2011                                                   #
# Version: 1.0                                                          #
#-----------------------------------------------------------------------#

## Main function
huge.NPN = function(x, NPN.func = "shrinkage", NPN.thresh = NULL, verbose = TRUE){
	gcinfo(FALSE)
	n = nrow(x)
  	d = ncol(x)
  	xt = list()
  	xt$NPN.func = NPN.func
  	xt$ntdata = x
  	
  	# Shrinkaage transformation
	if(NPN.func == "shrinkage"){
		if(verbose) cat("Conducting NonparaNormal (NPN) transformation via shrunkun ECDF....")
		
		x = qnorm(apply(x,2,rank)/(n+1))
		xt$data = x/sd(x[,1])
		
		if(verbose) cat("done.\n")
		rm(x,n,d,verbose)
   		gc()	
	}
	
	# Truncation transformation
	if(NPN.func == "truncation"){
		if(verbose) cat("Conducting NonparaNormal (NPN) transformation via truncated ECDF....")
		if(is.null(NPN.thresh)) NPN.thresh = 1/(4*(n^0.25)*sqrt(pi*log(n)))
		
		x = qnorm(pmin(pmax(apply(x,2,rank)/n, NPN.thresh), 1-NPN.thresh))
    	xt$data = x/sd(x[,1])
    	
    	if(verbose) cat("done.\n")
    	rm(x,n,d,NPN.thresh,verbose)
   		gc()
	}
	
	# Output class "NPN"
	class(xt) = "NPN"
	return(xt)
}

## Default print function for class "NPN" 
print.NPN = function(x, ...){
	cat("Gaussianized data by huge.NPN()\n")
	cat("Sample size: n =", nrow(x$data), "\n")
	cat("Dimension: d =", ncol(x$data), "\n")
	cat("NonparanNormal transformation type:", x$NPN.func,"\n")
}

## Default summary function for class "NPN"
summary.NPN = function(object, ...){
	cat("Gaussianized data by huge.NPN()\n")
	cat("Sample size: n =", nrow(object$data), "\n")
	cat("Dimension: d =", ncol(object$data), "\n")
	cat("NonparaNormal transformation type:", object$NPN.func,"\n")
}

## Default plot function for class "NPN"
plot.NPN = function(x, ...){
	par = par(mfrow = c(2,1), pty = "s", omi=c(0.5,0.5,0.5,0.5), mai = c(0.5,0.5,0.5,0.5))
	image(cor(x$ntdata), col = gray.colors(256), main = "Original Empirical Covariance Matrix")
	image(cor(x$data), col = gray.colors(256), main = "Transfored Empirical Covariance Matrix")
}