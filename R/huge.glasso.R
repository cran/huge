#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# glasso(): A Graphical Lasso (GLASSO) using Sparse Matrices           #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@gmail.com>; <hanliu@cs.jhu.edu>                     #
# Date: Jun 15th 2011                                                    #
# Version: 1.0.3                                                        #
#-----------------------------------------------------------------------#

## Main function
huge.glasso = function(x, lambda = NULL, lambda.min.ratio = NULL, nlambda = NULL, cov.output = FALSE, verbose = TRUE){
	
	gcinfo(FALSE)
	n = nrow(x)
	d = ncol(x)
	fit = list()
	fit$cov.input = isSymmetric(x)
	if(fit$cov.input)
	{
		if(verbose) cat("The input is identified as the covriance matrix.\n")
		S = x
	}
	if(!fit$cov.input)
	{
		x = scale(x)
		S = t(x)%*%x/n
	}
	rm(x)
	gc()
	
	if(!is.null(lambda)) nlambda = length(lambda)
	if(is.null(lambda))
	{
		if(is.null(nlambda))
			nlambda = 10
		if(is.null(lambda.min.ratio))
			lambda.min.ratio = 0.1
		lambda.max = max(max(S-diag(d)),-min(S-diag(d)))
		lambda.min = lambda.min.ratio*lambda.max
		lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
	}
	
	fit$lambda = lambda
	fit$loglik = rep(0,nlambda)
	fit$sparsity = rep(0,nlambda)
	fit$df = rep(0,nlambda)
	if(verbose)
	{
		cat(paste(c("Conducting Graphical Lasso (GLASSO)....in progress:", floor(100*((1/nlambda)^2)), "%"), collapse=""), "\r")
		flush.console()
	}
	
	out.glasso = glasso(S, lambda[nlambda])
	tmp.w = out.glasso$w
	tmp.wi = out.glasso$wi
	fit$loglik[nlambda] = out.glasso$loglik + lambda[nlambda]*sum(abs(tmp.wi[[nlambda]]))
	rm(out.glasso)
	gc()
	fit$wi = list()
	fit$wi[[nlambda]] = Matrix(tmp.wi,sparse = TRUE)
	if(cov.output)
	{
		fit$w = list()
		fit$w[[nlambda]] = Matrix(tmp.w,sparse = TRUE)
	}
	
	fit$path = list()
	fit$path[[nlambda]] = Matrix(abs(sign(tmp.wi)), sparse = TRUE)
	diag(fit$path[[nlambda]]) = 0
	
	fit$df[nlambda] = sum(fit$path[[nlambda]])/2
	fit$sparsity[nlambda] = 2*fit$df[[nlambda]]/d/(d-1)
	if(nlambda>1)
		for(i in (nlambda-1):1){
			if(verbose){
				cat(paste(c("Conducting Graphical Lasso (GLASSO)....in progress:", floor(100*(((nlambda-i)/nlambda)^2)), "%"), collapse=""), "\r")
				flush.console()
			}
			out.glasso = glasso(S,rho = fit$lambda[i], w.init = tmp.w, wi.init = tmp.wi, start = "warm")
			tmp.w = out.glasso$w
			tmp.wi = out.glasso$wi
			fit$loglik[i] = out.glasso$loglik + lambda[i]*sum(abs(tmp.wi))
			rm(out.glasso)
			gc()
			fit$wi[[i]] = Matrix(tmp.wi,sparse = TRUE)
			if(cov.output)
				fit$w[[i]] = Matrix(tmp.w,sparse = TRUE)
			fit$path[[i]] = Matrix(abs(sign(tmp.wi)), sparse = TRUE)
			diag(fit$path[[i]]) = 0
			fit$df[i] = sum(fit$path[[i]])/2
			fit$sparsity[i] = 2*fit$df[[i]]/d/(d-1)
			gc()
		}
	rm(S,tmp.w,tmp.wi)
	gc()
	if(verbose){
   		cat("Conducting Graphical Lasso (GLASSO)....done.              \r")
   		cat("\n")
        flush.console()
    }
    return(fit)
}