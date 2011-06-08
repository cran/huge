#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.GECT():                                                          #                           
# Graph Estimation via Correlation Thresholding (GECT)                  # 
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@gmail.com>; <hanliu@cs.jhu.edu>                     #
# Date: Jun 8th 2011                                                    #
# Version: 1.0.2                                                        #
#-----------------------------------------------------------------------#

##Main function
huge.GECT = function(x, nlambda = NULL, lambda.min.ratio = NULL, lambda = NULL, verbose = TRUE)
{	
	gcinfo(FALSE)
  	d = ncol(x)
  	fit = list()
  	fit$cov.input = isSymmetric(x)
  	if(!fit$cov.input)	S = abs(cor(x))
  	if(fit$cov.input)
  	{
  		if(verbose) cat("The input is identified as the covriance matrix.\n")
  		S = abs(x);
  	}
  	diag(S) = 0
  	S.rank = order(S,decreasing = TRUE)
  	rm(x)
	gc()
 		
	if(is.null(lambda))
	{
		if(is.null(nlambda))
			nlambda = 20
		if(is.null(lambda.min.ratio))
			lambda.min.ratio = 0.05
		 		
 		density.max = lambda.min.ratio*d*(d-1)/2
 		density.min = 1
 		density.all = ceiling(seq(density.min,density.max,length = nlambda))*2
 		fit$sparsity = density.all/d/(d-1)
 		fit$lambda = S[S.rank[density.all]]
 		rm(density.max,lambda.min.ratio,density.min,S)
		gc()
 						
 		fit$path = list()
		for(i in 1:nlambda)
		{
			fit$path[[i]] = Matrix(0,d,d)
			fit$path[[i]][S.rank[1:density.all[i]]] = 1
			if(verbose)
			{
   				cat(paste(c("Conducting Graph Estimation via Correlation Thresholding (GECT)....in progress:", floor(100*i/nlambda), "%"), collapse=""), "\r")
            	flush.console()
            }	
		}
		rm(density.all,nlambda,S.rank)
		gc()
	}

	if(!is.null(lambda))
	{
		nlambda = length(lambda)
		fit$path = list()
		fit$sparsity = rep(0,nlambda)
		for(i in 1:nlambda)
		{
			fit$path[[i]] = Matrix(0,d,d)
			fit$path[[i]][S > lambda[i]] = 1
			fit$sparsity[i] = sum(fit$path[[i]])/d/(d-1)
			if(verbose)
			{
   				mes <- paste(c("Conducting Graph Estimation via Correlation Thresholding (GECT)....in progress:", floor(100*i/nlambda), "%"), collapse="")
   				cat(mes, "\r")
            	flush.console()
            }
		}
		fit$lambda = lambda
		rm(S,lambda)
		gc()
	}
		
	if(verbose)
	{
        cat("Conducting Graph Estimation via Correlation Thresholding (GECT)....done.             \r\n")
        flush.console()
    }
	gc()
	class(fit) = "GECT"
	return(fit)
}

# Default printing function
print.GECT = function(x, ...)
{
	cat("This is a solution path using Graph Estimation via Correlation Thresholding (GECT) and length = ", length(x$path), "\n")
	cat("huge.GECT() is an internal function, please refer to huge() and huge.select()")
}

# Default plot function
plot.GECT = function(x, ...)
{
	par(mfrow = c(1,1),pty = "s")
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "b",xlim = rev(range(x$lambda)))
}