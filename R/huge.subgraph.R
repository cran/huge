#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.subgraph(): M&B Graph Estimation via Lasso (GEL)                 #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 12th 2010                                                   #
# Version: 0.8.1                                                         #
#-----------------------------------------------------------------------#

## Main function
huge.subgraph = function(x, ind.group = NULL, ind.mat = NULL, alpha = 1, lambda = NULL, nlambda = 10, lambda.min.ratio = 0.1, sym = "or", verbose = TRUE){

	gcinfo(FALSE)
	fit = list()
	n = nrow(x)
	d = ncol(x)
	
	if(d < 3){
		cat("The fullgraph dimension < 3 and huge.subgraph() will be teminated....\n")
		cat("Please refer to Pearson's product-moment correlation....\n")
		fit$marker = "Terminated"
		class(fit) = "subgraph"
		return(fit)
	}
		
	if(!is.null(ind.mat))
		if(!is.matrix(ind.mat)){
			cat("The neighborhood size < 2 and huge.subgraph() will be teminated.\n")
			cat("Please refer to Pearson's product-moment correlation....\n")
			fit$marker = "Terminated"
			class(fit) = "subgraph"
			return(fit)
	}
	
	if(is.null(ind.group))	ind.group = c(1:d)
  	k = length(ind.group)
		
	x = scale(x)
	if(is.null(lambda)){
		if(!is.null(ind.mat)){
			tmp = rep(0,k)
			for (i in 1:k)	tmp[i] = max(abs(t(x[,ind.mat[,i]])%*%x[,ind.group[i]]))
			lambda = max(tmp)/n/alpha*exp(seq(log(1), log(lambda.min.ratio), length = nlambda))
			rm(tmp)
			gc()
		}
		if(is.null(ind.mat)){
			tmp = rep(0,k)
			for (i in 1:k)	tmp[i] = max(abs(t(x[,-ind.group[i]])%*%x[,ind.group[i]]))
			lambda = max(tmp)/n/alpha*exp(seq(log(1), log(lambda.min.ratio), length = nlambda))
			rm(tmp)
			gc()
		}
	}
	fit$lambda = lambda
	nlambda = length(lambda)

	fit$rss = matrix(0,k,nlambda)
   	fit$df = matrix(0,k,nlambda)

   	fit$path = list()
   	for(i in 1:nlambda)	fit$path[[i]] = Matrix(0,k,k)
   	if(is.null(ind.mat)){
   		
   		for(i in 1:k){
   			if(verbose&&alpha == 1)
   				mes <- paste(c("Conducting Meinshausen & Buhlmann Graph Estimation via Lasso (GEL)....in progress:", floor(100*i/k), "%"), collapse="")
   			if(verbose&&alpha < 1)
   				mes <- paste(c("Conducting Meinshausen & Buhlmann Graph Estimation via Elastic Net....in progress:", floor(100*i/k), "%"), collapse="")
   			if(verbose){
   				cat(mes, "\r")
            	flush.console()
            }
   			out.glm = glmnet(x[,-ind.group[i]], x[,ind.group[i]], lambda = lambda, alpha = alpha)
      		fit$rss[i,] = out.glm$nulldev*(1-out.glm$dev)
      		fit$df[i,] = out.glm$df
     		
      		for(j in 1:nlambda){
         		tmp = rep(0,d)
    	 		tmp[-ind.group[i]] = sign(abs(out.glm$beta[,j])!=0)
    	  		fit$path[[j]][i,] = tmp[ind.group]
      		}
      		#Sys.sleep(0.01)
   			#if(verbose)	setTxtProgressBar(pb, i)
      	}
   		rm(x,n,d,ind.group,lambda,tmp,out.glm)
   		gc()
   	}
   	
   	if(!is.null(ind.mat)){
   		for(i in 1:k){
   			if(verbose&&alpha == 1)
   				mes <- paste(c("Conducting Meinshausen & Buhlmann Graph Estimation via Lasso (GEL)....in progress:", floor(100*i/k), "%"), collapse="")
   			if(verbose&&alpha < 1)
   				mes <- paste(c("Conducting Meinshausen & Buhlmann Graph Estimation via Elastic Net....in progress:", floor(100*i/k), "%"), collapse="")
   			if(verbose){
   				cat(mes, "\r")
            	flush.console()
            }
      		out.glm = glmnet(x[,ind.mat[,i]], x[,ind.group[i]], lambda = lambda, alpha = alpha)
      		
      		fit$rss[i,] = out.glm$nulldev*(1-out.glm$dev)
      		fit$df[i,] = out.glm$df
      		
      		for(j in 1:nlambda){
         		tmp = rep(0,d)
    	 		tmp[ind.mat[,i]] = sign(abs(out.glm$beta[,j])!=0)
    	  		fit$path[[j]][i,] = tmp[ind.group]
      		}
      		#Sys.sleep(0.01)
   			#if(verbose)	setTxtProgressBar(pb, i)
   		}
   		rm(x,n,d,ind.group,ind.mat,lambda,tmp,out.glm)
   		gc()
   	}
   	#if(verbose)	close(pb)
	if(verbose&&alpha == 1)
   		mes <- "Conducting Meinshausen & Buhlmann Graph Estimation via Lasso (GEL)....done.              "
   	if(verbose&&alpha < 1)
   		mes <- "Conducting Meinshausen & Buhlmann Graph Estimation via Elastic Net....done.              "
   	if(verbose){
   		cat(mes, "\r")
   		cat("\n")
        flush.console()
    }
	fit$sparsity = rep(0,nlambda)
   	if(sym == "or")
   		for(i in 1:nlambda){
   			fit$path[[i]]= sign(fit$path[[i]] + t(fit$path[[i]]))
   			fit$sparsity[i] = sum(fit$path[[i]])/k/(k-1)
   		}
   	if(sym == "and")
   		for(i in 1:nlambda){
   			fit$path[[i]]= fit$path[[i]]*t(fit$path[[i]])
   			fit$sparsity[i] = sum(fit$path[[i]])/k/(k-1)
   		}
   		
   	rm(verbose,nlambda,k)
   	gc()
   	fit$marker = "Successful"
   	class(fit) = "subgraph"
   	return(fit)
}

# Default printing function
print.subgraph = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.subgraph() has been terminated\n")
		return("Please refer to the manual")
	}
	cat("This is a solution path using Meinshausen & Buhlman Graph Estimation via Lasso or Elastic Net and length = ", length(x$path), "\n")
	cat("huge.subgraph() is an internal function. For more information, please refer to huge() and huge.select().\n")
}
	
# Defaulty summary function	
summary.subgraph = function(object, ...){
	if(object$marker == "Terminated"){
		cat("huge.subgraph() has been terminated\n")
		return("Please refer to the manual")
	}
	cat("This is a solution path using Meinshausen & Buhlman Graph Estimation via Lasso or Elastic Net and length = ", length(object$path), "\n")
	cat("huge.subgraph() is an internal function. For more information, please refer to huge() and huge.select().\n")
}

# Default plot function
plot.subgraph = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.subgraph() has been terminated\n")
		return("Please refer to the manual")
	}
	par(mfrow = c(1,1))
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameters", ylab = "Sparsity Level", type = "b",xlim = rev(range(x$lambda)))
	cat("huge.subgraph() is an internal function. For more information, please refer to huge() and huge.select().\n")
}