#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.npn(): NonparaNormal transofmration                              #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 12th, 2010                                                  #
# Version: 0.8                                                          #
#-----------------------------------------------------------------------#

## Main function
huge.npn = function(x, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE){
	
	n = nrow(x)
  	d = ncol(x)
  	xt = list()
  	xt$npn.func = npn.func
  	xt$ntdata = x
  	
  	# Shrinkaage transformation
	if(npn.func == "shrinkage"){
		if(verbose) cat("Conducting NonparaNormal (NPN) transformation via shrunkun ECDF....")
		for (i in 1:d) x[,i] = qnorm(rank(x[,i])/(n + 1))
		xt$data = x/sd(x[,1])
		if(verbose) cat("done.\n")
		rm(n, d, verbose)
   		gc(gcinfo(verbose = FALSE)) 	
	}
	
	# Truncation transformation
	if(npn.func == "truncation"){
		if(verbose) cat("Conducting NonparaNormal (NPN) transformation via truncated ECDF....")
		if(is.null(npn.thresh)) npn.thresh = 1/(4*(n^0.25)*sqrt(pi*log(n)))
		for (i in 1:d) {
    		tmp = rank(x[,i])/n
    		tmp = pmax(tmp, npn.thresh)
    		tmp = pmin(tmp, 1-npn.thresh)
       		x[,i] = qnorm(tmp)
    	}
    	xt$data = x/sd(x[,1])
    	
    	if(verbose) cat("done.\n")
    	
    	rm(n,d,tmp,npn.func,npn.thresh,verbose)
   		gc(gcinfo(verbose = FALSE))
	}
	
	# Output class "npn"
	class(xt) = "npn"
	return(xt)
}

## Default print function for class "npn" 
print.npn = function(x, ...){
	cat("Gaussianized data by huge.npn()\n")
	cat("Sample size: n =", nrow(x$data), "\n")
	cat("Dimension: d =", ncol(x$data), "\n")
	cat("NonparanNormal transformation type:", x$npn.func,"\n")
}

## Default summary function for class "npn"
summary.npn = function(object, ...){
	cat("Gaussianized data by huge.npn()\n")
	cat("Sample size: n =", nrow(object$data), "\n")
	cat("Dimension: d =", ncol(object$data), "\n")
	cat("NonparaNormal transformation type:", object$npn.func,"\n")
}

## Default plot function for class "npn"
plot.npn = function(x, ...){
	par = par(mfrow = c(2,1), pty = "s", omi=c(0.5,0.5,0.5,0.5), mai = c(0.5,0.5,0.5,0.5))
	image(cor(x$ntdata), col = gray.colors(256), main = "Original Empirical Covariance Matrix")
	image(cor(x$data), col = gray.colors(256), main = "Transfored Empirical Covariance Matrix")
}