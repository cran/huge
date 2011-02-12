#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# glassoM(): A sligthly modifeid Graphical Lasso (GLASSO)              #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Feb 12th 2011                                                   #
# Version: 0.9.1                                                          #
#-----------------------------------------------------------------------#

## Main function
huge.glassoM = function(x, ind.group = NULL, lambda.min.ratio = 0.1, nlambda = 10, lambda = NULL, verbose = TRUE){
	
	gcinfo(FALSE)
	n = nrow(x)
	d = ncol(x)
	fit = list()
	
	if(d < 3){
		cat("The fullgraph dimension < 3 and huge.subgraph() will be teminated....\n")
		cat("Please refer to Pearson's product-moment correlation....\n")
		rm(x,lambda.min.ratio,nlambda)
		gc()
		fit$marker = "Terminated"
		class(fit) = "glassoM"
		return(fit)
	}
	
	s = cor(x)
	rm(x)
	gc()
	
	if(is.null(ind.group)) ind.group = 1:d
	k = length(ind.group)
	tmp = max(abs(s-diag(d)))
	if(is.null(lambda)){
		lambda = exp(seq(log(tmp),log(lambda.min.ratio*tmp),length = nlambda))
	}
	rm(tmp)
	gc()
	nlambda = length(lambda)
	fit$lambda = lambda
	fit$loglik = rep(0,nlambda)
	fit$sparsity = rep(0,nlambda)
	fit$df = rep(0,nlambda)
	if(verbose){
		cat(paste(c("Conducting Graphical Lasso (GLASSO)....in progress:", floor(100*1/nlambda), "%"), collapse=""), "\r")
		flush.console()
	}
	out.glasso = glasso(s, lambda[1])
	fit$wi = list()
	fit$wi[[1]] = Matrix(out.glasso$wi,sparse = TRUE)
	tmp.w = out.glasso$w
	diag(out.glasso$wi) = 0
	fit$path = list()	
	fit$path[[1]] = Matrix(abs(sign(out.glasso$wi[ind.group,ind.group])), sparse = TRUE)
	fit$loglik[1] = out.glasso$loglik + lambda[1]*sum(abs(fit$wi[[1]]))
	fit$df[1] = sum(abs(sign(out.glasso$wi)))/2
	fit$sparsity[1] = 2*fit$df[[1]]/k/(k-1)
	rm(out.glasso)
	gc()
	if(nlambda>1)
		for(i in 2:nlambda){
			if(verbose){
				cat(paste(c("Conducting Graphical Lasso (GLASSO)....in progress:", floor(100*((i/nlambda)^2)), "%"), collapse=""), "\r")
				flush.console()
			}
			out.glasso = glasso(s,rho = fit$lambda[i], start = "warm", w.init = tmp.w, wi.init = as.matrix(fit$wi[[i-1]]))
			tmp.w = out.glasso$w
			fit$wi[[i]] = Matrix(out.glasso$wi,sparse = TRUE)
			diag(out.glasso$wi) = 0
			fit$path[[i]] = Matrix(abs(sign(out.glasso$wi[ind.group,ind.group])), sparse = TRUE)
			fit$loglik[i] = out.glasso$loglik + lambda[i]* sum(abs(fit$wi[[i]]))
			fit$df[i] = sum(abs(sign(out.glasso$wi)))/2
			fit$sparsity[i] = sum(fit$path[[i]])/k/(k-1)
			rm(out.glasso)
			gc()
		}
	rm(s,tmp.w)
	gc()
	if(verbose){
   		cat("Conducting Graphical Lasso (GLASSO)....done.              \r")
   		cat("\n")
        flush.console()
    }
    fit$marker = "Successful"
    class(fit) = "glassoM"
    return(fit)
}

# Default printing function
print.glassoM = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.glassoM() has been terminated\n")
		return("Please refer to the manual")
	}
	cat("This is a solution path using Graphical Lasso (GLASSO) and length = ", length(x$path), "\n")
	cat("huge.glassoM() is an internal function. For more information, please refer to huge() and huge.select().\n")
}
	
# Defaulty summary function	
summary.glassoM = function(object, ...){
	if(object$marker == "Terminated"){
		cat("huge.glassoM() has been terminated\n")
		return("Please refer to the manual")
	}
	cat("This is a solution path using Graphical Lasso (GLASSO) and length = ", length(object$path), "\n")
	cat("huge.glassoM() is an internal function. For more information, please refer to huge() and huge.select().\n")
}

# Default plot function
plot.glassoM = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.glassoM() has been terminated\n")
		return("Please refer to the manual")
	}
	par(mfrow = c(1,1))
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "b",xlim = rev(range(x$lambda)))
	cat("huge.glassoM() is an internal function. For more information, please refer to huge() and huge.select().\n")
}