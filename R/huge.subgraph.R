#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.subgraph(): subgraph estimation                                  #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@gmail.com>; <hanliu@cs.jhu.edu>                     #
# Date: Nov 9th 2010                                                    #
# Version: 0.7                                                          #
#-----------------------------------------------------------------------#
# Main function
huge.subgraph = function(x, ind.group = NULL, ind.mat = NULL, alpha = 1, lambda = NULL, n.lambda = 10, lambda.min = 0.1, sym = "or", verbose = TRUE){

	n = nrow(x)
	d = ncol(x)
	fit = list()
	if(is.null(ind.group))	ind.group = c(1:d)
  	k = length(ind.group)
  	
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
		
	x = scale(x)
	if(is.null(lambda)){
		if(!is.null(ind.mat)){
			tmp = rep(0,k)
			for (i in 1:k)	tmp[i] = max(abs(t(x[,ind.mat[,i]])%*%x[,ind.group[i]]))
			lambda = max(tmp)/n/alpha*exp(seq(log(1), log(lambda.min), length = n.lambda))
			rm(tmp)
			gc(gcinfo(verbose = FALSE))
		}
		if(is.null(ind.mat)){
			tmp = rep(0,k)
			for (i in 1:k)	tmp[i] = max(abs(t(x[,-ind.group[i]])%*%x[,ind.group[i]]))
			lambda = max(tmp)/n/alpha*exp(seq(log(1), log(lambda.min), length = n.lambda))
			rm(tmp)
			gc(gcinfo(verbose = FALSE))
		}
	}
	fit$lambda = lambda
	n.lambda = length(lambda)

	fit$rss = matrix(0,k,n.lambda)
   	fit$df = matrix(0,k,n.lambda)
   	fit$sigma2hat = rep(0,k)

   	fit$path = list()
   	for(i in 1:n.lambda)	fit$path[[i]] = Matrix(0,k,k)
    
    if(verbose) cat("Computing the solution path:")
   	pdot = ceiling(k/30)
   	if(is.null(ind.mat)){
   		for(i in 1:k){
   			if(verbose&&(i%%pdot) == 1) cat(".")
      		out.glm = glmnet(x[,-ind.group[i]], x[,ind.group[i]], lambda = lambda, alpha = alpha)
      		
      		fit$rss[i,] = out.glm$nulldev*(1-out.glm$dev)
      		fit$df[i,] = out.glm$df
     		
      		for(j in 1:n.lambda){
         		tmp = rep(0,d)
    	 		tmp[-ind.group[i]] = sign(abs(out.glm$beta[,j])!=0)
    	  		fit$path[[j]][i,] = tmp[ind.group]
      		}
      	}
   		rm(x,n,d,ind.group,alpha,lambda,tmp,out.glm)
   		gc(gcinfo(verbose = FALSE))
   	}
   	if(!is.null(ind.mat)){
   		for(i in 1:k){
   			if(verbose&&(i%%pdot) == 1) cat(".")
      		out.glm = glmnet(x[,ind.mat[,i]], x[,ind.group[i]], lambda = lambda, alpha = alpha)
      		
      		fit$rss[i,] = out.glm$nulldev*(1-out.glm$dev)
      		fit$df[i,] = out.glm$df
      		
      		for(j in 1:n.lambda){
         		tmp = rep(0,d)
    	 		tmp[ind.mat[,i]] = sign(abs(out.glm$beta[,j])!=0)
    	  		fit$path[[j]][i,] = tmp[ind.group]
      		}
   		}
   		rm(x,n,d,ind.group,ind.mat,alpha,lambda,tmp,out.glm)
   		gc(gcinfo(verbose = FALSE))
   	}
   	if(verbose)	cat("done.\n")
	
	fit$sparsity = rep(0,n.lambda)
   	if(sym == "or")
   		for(i in 1:n.lambda){
   			fit$path[[i]]= sign(fit$path[[i]] + t(fit$path[[i]]))
   			fit$sparsity[i] = sum(fit$path[[i]])/k/(k-1)
   		}
   	if(sym == "and")
   		for(i in 1:n.lambda){
   			fit$path[[i]]= fit$path[[i]]*t(fit$path[[i]])
   			fit$sparsity[i] = sum(fit$path[[i]])/k/(k-1)
   		}
   		
   	rm(verbose,n.lambda,k)
   	gc(gcinfo(verbose = FALSE))
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
	cat("This is a solution path using Meinshausen Buhlman Graph Estimation via Lasso or elastic net and length = ", length(x$path), "\n")
	cat("huge.subgraph is an internal function in the package, please refer to huge() and huge.select().\n")
}
	
# Defaulty summary function	
summary.subgraph = function(object, ...){
	if(object$marker == "Terminated"){
		cat("huge.subgraph() has been terminated\n")
		return("Please refer to the manual")
	}
	cat("This is a solution path using Meinshausen & Buhlman Graph Estimation via Lasso or elastic net and length = ", length(object$path), "\n")
	cat("huge.subgraph is an internal function in the package, please refer to huge() and huge.select().\n")
}

# Default plot function
plot.subgraph = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.subgraph() has been terminated\n")
		return("Please refer to the manual")
	}
	par(mfrow = c(1,1))
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameters", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)))
	cat("huge.subgraph is an internal function used in the package, please refer to huge() and huge.select().\n")
}