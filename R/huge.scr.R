#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.scr():                                                           #
# (1) Graph SURE Screening (GSS)                                        #
# (2) Graph Approximation via Correlation Thresholding (GACT)           # 
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 12th 2010                                                   #
# Version: 0.8                                                          #
#-----------------------------------------------------------------------#

##Main function
huge.scr = function(x, ind.group = NULL, scr.num = NULL, approx = FALSE, nlambda = 30, lambda.min.ratio = 0.05, lambda = NULL, verbose = TRUE){
	
	n = nrow(x)
  	d = ncol(x)
  	fit = list()
  	
  	if(is.null(ind.group))	ind.group = c(1:d)
  	k = length(ind.group)
  	
  	if(d < 3){
		cat("The fullgraph dimension < 3 and huge.scr() will be teminated....\n")
		cat("Please refer to Pearson's product-moment correlation....\n")
		fit$marker = "Terminated"
		class(fit) = "scr"
		return(fit)
	}
	
	if(!is.null(scr.num))
		if(scr.num == 1){
			cat("The neighborhood size < 2 and huge.scr() will be teminated.\n")
			cat("Please refer to Pearson's product-moment correlation....\n")
			fit$marker = "Terminated"
			class(fit) = "scr"
			return(fit)
		}
  	
  	fit$approx = approx
  	
  	
  	if(approx){
  		
  		xt = x[,ind.group]
  		rm(x)
		gc(gcinfo(verbose = FALSE))	
		xt = scale(xt)
 				  		
  		S = abs(cov(xt))
  		diag(S) = 0
		lambda.max = max(S) 
		lambda.min = lambda.min.ratio*lambda.max
		
		lambda = 1-exp(seq(log(1-lambda.max),log(1-lambda.min),length = nlambda))
		rm(lambda.max,lambda.min.ratio,lambda.min)
		gc(gcinfo(verbose = FALSE))
			
		fit$path = list()
		for(i in 1:nlambda)	fit$path[[i]] = Matrix(0,k,k)
		for(i in 1:nlambda){
			if(verbose){
   				mes <- paste(c("Conducting Graph Approximation via Correlation Thresholding (GACT)....in progress:", floor(100*i/nlambda), "%"), collapse="")
   				cat(mes, "\r")
            	flush.console()
            }
			fit$path[[i]][S > lambda[i]] = 1
		}
		if(verbose){
			mes = "Conducting Graph Approximation via Correlation Thresholding (GACT)....done.             "
        	cat(mes, "\r")
        	cat("\n")
        	flush.console()
        }
		
		fit$sparsity = rep(0,nlambda)
		for(i in 1:nlambda) fit$sparsity[i] = sum(fit$path[[i]])/k/(k-1)
		
		fit$lambda = lambda
		rm(lambda,xt)
		gc(gcinfo(verbose = FALSE))

	}
	
	if(!approx){
		if(is.null(scr.num))	scr.num = min(n-1,d-1)
		if(verbose) cat("Conducting Graph SURE Screening (GSS)....")
		
		fit$ind.mat = matrix(0,scr.num,k)
		
		# Patition ind.group
		g = ceiling(d*k/1e6)
		l = floor(k/g)
		r = k%%l
	
		# Compute correlation blocks
		if(g>0)
			for(i in 1:g){
				cor.block = -abs(t(x)%*%x[,ind.group[((i-1)*l+1):(i*l)]])
				fit$ind.mat[,((i-1)*l+1):(i*l)] = apply(cor.block,2,order)[2:(scr.num+1),]
				rm(cor.block)
				gc(gcinfo(verbose = FALSE))
			}
		
		if(r>0){
			cor.block = -abs(t(x)%*%x[,ind.group[(g*l+1):(g*l+r)]])
			fit$ind.mat[,(g*l+1):(g*l+r)] = apply(cor.block,2,order)[2:(scr.num+1),]
			rm(cor.block)
		}
		rm(x)
		if(verbose) cat("done.\n")
	}
	
	fit$marker = "Successful"
	
	gc(gcinfo(verbose = FALSE))
	class(fit) = "scr"
	return(fit)
}

# Default printing function
print.scr = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.scr() has been terminated\n")
		return("Please refer to the manual")
	}
	if(x$approx) cat("This is a solution path using Graph Approximation via Correlation Thresholding (GACT) and length = ", length(x$path), "\n")
	if(!x$approx) cat("The dimension of prelected neighborhood after screening = ", nrow(x$ind.mat), "\n")
	cat("huge.scr() is an internal function, please refer to huge() and huge.select()")
}

# Default summary function
summary.scr = function(object, ...){
	if(object$marker == "Terminated"){
		cat("huge.scr() has been terminated\n")
		return("Please refer to the manual")
	}
	if(object$approx) cat("This is a solution path using Graph Approximation via Correlation Thresholding (GACT) and length = ", length(object$path), "\n")
	if(!object$approx) cat("The dimension of preselected neighborhood after screening = ", nrow(object$ind.mat), "\n")
	cat("huge.scr() is an internal function, please refer to huge() and huge.select()")
}

# Default plot function
plot.scr = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.scr() has been terminated\n")
		return("Please refer to the manual")
	}
	if(x$approx){
		par(mfrow = c(1,1))
		tmp = (x$lambda > 0)
		plot(x$lambda[tmp], x$sparsity[tmp], log = "x", xlab = "Regularization Parameters", ylab = "Sparsity Level", type = "b",xlim = rev(range(x$lambda)))
	}
	if(!x$approx) cat("No plot information avaiable for the Graph SURE Screening (GSS)","\n")
	cat("huge.scr() is an internal function. For more information, please refer to huge() and huge.select()")
}