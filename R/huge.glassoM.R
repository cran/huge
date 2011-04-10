#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# glassoM(): A Graphical Lasso (GLASSO) using Sparse Matrices           #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Apr 10th 2011                                                   #
# Version: 1.0.1                                                        #
#-----------------------------------------------------------------------#

## Main function
huge.glassoM = function(x, lambda = NULL, lambda.min.ratio = NULL, nlambda = NULL, cov.glasso = FALSE, verbose = TRUE){
	
	gcinfo(FALSE)
	n = nrow(x)
	d = ncol(x)
	fit = list()
	
	S = t(x)%*%x/n;
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
		cat(paste(c("Conducting Graphical Lasso (GLASSO)....in progress:", floor(100*1/nlambda), "%"), collapse=""), "\r")
		flush.console()
	}
	
	out.glasso = glasso(S, lambda[nlambda])
	fit$wi = list()
	fit$wi[[nlambda]] = Matrix(out.glasso$wi,sparse = TRUE)
	if(cov.glasso)
	{
		fit$w = list()
		fit$w[[nlambda]] = Matrix(out.glasso$w,sparse = TRUE)
	}
	
	tmp.w = out.glasso$w
		
	diag(out.glasso$wi) = 0
	fit$path = list()
	fit$path[[nlambda]] = Matrix(abs(sign(out.glasso$wi)), sparse = TRUE)
	fit$loglik[nlambda] = out.glasso$loglik + lambda[nlambda]*sum(abs(fit$wi[[nlambda]]))
	fit$df[nlambda] = sum(fit$path[[nlambda]])/2
	fit$sparsity[nlambda] = 2*fit$df[[nlambda]]/d/(d-1)
	rm(out.glasso)
	gc()
	if(nlambda>1)
		for(i in (nlambda-1):1){
			if(verbose){
				cat(paste(c("Conducting Graphical Lasso (GLASSO)....in progress:", floor(100*(((nlambda-i)/nlambda)^2)), "%"), collapse=""), "\r")
				flush.console()
			}
			out.glasso = glasso(S,rho = fit$lambda[i], w.init = tmp.w, wi.init = as.matrix(fit$wi[[i+1]]))
			tmp.w = out.glasso$w
			fit$wi[[i]] = Matrix(out.glasso$wi,sparse = TRUE)
			if(cov.glasso)
				fit$w[[i]] = Matrix(out.glasso$w,sparse = TRUE)
			diag(out.glasso$wi) = 0
			fit$path[[i]] = Matrix(abs(sign(out.glasso$wi)), sparse = TRUE)
			fit$loglik[i] = out.glasso$loglik + lambda[i]* sum(abs(fit$wi[[i]]))
			fit$df[i] = sum(abs(sign(out.glasso$wi)))/2
			fit$sparsity[i] = sum(fit$path[[i]])/d/(d-1)
			rm(out.glasso)
			gc()
		}
	rm(S,tmp.w)
	gc()
	if(verbose){
   		cat("Conducting Graphical Lasso (GLASSO)....done.              \r")
   		cat("\n")
        flush.console()
    }
    class(fit) = "glassoM"
    return(fit)
}

# Default printing function
print.glassoM = function(x, ...){
	cat("This is a solution path using Graphical Lasso (GLASSO) and length = ", length(x$path), "\n")
	cat("huge.glassoM() is an internal function. For more information, please refer to huge() and huge.select().\n")
}
	
# Defaulty summary function	
summary.glassoM = function(object, ...){
	cat("This is a solution path using Graphical Lasso (GLASSO) and length = ", length(object$path), "\n")
	cat("huge.glassoM() is an internal function. For more information, please refer to huge() and huge.select().\n")
}

# Default plot function
plot.glassoM = function(x, ...){
	par(mfrow = c(1,1))
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "b",xlim = rev(range(x$lambda)))
	cat("huge.glassoM() is an internal function. For more information, please refer to huge() and huge.select().\n")
}