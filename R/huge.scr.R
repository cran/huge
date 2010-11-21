#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.scr():                                                           #
# (1) Graph SURE Screening (GSS)                                        #
# (2) Graph Approximation via Correlation Thresholding (GACT)           # 
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 21st 2010                                                   #
# Version: 0.9                                                          #
#-----------------------------------------------------------------------#

##Main function
huge.scr = function(x, ind.group = NULL, scr.num = NULL, method = "GSS", nlambda = 30, lambda.min.ratio = 0.1, lambda = NULL, verbose = TRUE){
	
	gcinfo(FALSE)
	n = nrow(x)
  	d = ncol(x)
  	fit = list()
  	
  	if(d < 3){
		cat("The fullgraph dimension < 3 and huge.scr() will be teminated....\n")
		cat("Please refer to Pearson's product-moment correlation....\n")
		rm(x,nlambda,lambda.min.ratio,n,d)
		gc()
		fit$marker = "Terminated"
		class(fit) = "scr"
		return(fit)
	}
	
	if(!is.null(scr.num))
		if(scr.num == 1){
			cat("The neighborhood size < 2 and huge.scr() will be teminated.\n")
			cat("Please refer to Pearson's product-moment correlation....\n")
			rm(x,nlambda,lambda.min.ratio,n,d)
			gc()
			fit$marker = "Terminated"
			class(fit) = "scr"
			return(fit)
		}
  	
  	if(is.null(ind.group))	ind.group = c(1:d)
  	k = length(ind.group)
  	
  	if(method == "GACT"){
  		S = abs(cor(x[,ind.group]))
  		diag(S) = 0
  		S.rank = order(S,decreasing = TRUE)
  		rm(x)
		gc()
 		
 		if(is.null(lambda)){
 			if(is.null(nlambda)) nlambda = 30
 			density.max = lambda.min.ratio*k*(k-1)/2
 			density.min = 1
 			density.all = ceiling(seq(density.min,density.max,length = nlambda))*2
 			fit$sparsity = density.all/k/(k-1)
 			fit$lambda = S[S.rank[density.all]]
 			rm(density.max,lambda.min.ratio,density.min,S)
			gc()
 			
 					
 			fit$path = list()
			for(i in 1:nlambda){
				fit$path[[i]] = Matrix(0,k,k)
				fit$path[[i]][S.rank[1:density.all[i]]] = 1
				if(verbose){
   					cat(paste(c("Conducting Graph Approximation via Correlation Thresholding (GACT)....in progress:", floor(100*i/nlambda), "%"), collapse=""), "\r")
            		flush.console()
            	}	
			}
			rm(density.all,nlambda,S.rank)
			gc()
		}

		if(!is.null(lambda)){
			nlambda = length(lambda)
			fit$path = list()
			fit$sparsity = rep(0,nlambda)
			for(i in 1:nlambda){
				fit$path[[i]] = Matrix(0,k,k)
				fit$path[[i]][S > lambda[i]] = 1
				fit$sparsity[i] = sum(fit$path[[i]])/k/(k-1)
				if(verbose){
   					mes <- paste(c("Conducting Graph Approximation via Correlation Thresholding (GACT)....in progress:", floor(100*i/nlambda), "%"), collapse="")
   					cat(mes, "\r")
            		flush.console()
            	}
			}
			fit$lambda = lambda
			rm(S,lambda)
			gc()
		}
		
		if(verbose){
        	cat("Conducting Graph Approximation via Correlation Thresholding (GACT)....done.             \r\n")
        	flush.console()
        }
		rm(n,d,k)
		gc()
	}
	
	if(method == "GSS"){
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
				gc()
			}
		
		if(r>0){
			cor.block = -abs(t(x)%*%x[,ind.group[(g*l+1):(g*l+r)]])
			fit$ind.mat[,(g*l+1):(g*l+r)] = apply(cor.block,2,order)[2:(scr.num+1),]
			rm(cor.block)
		}
		rm(x,r,g,l,scr.num,lambda.min.ratio,nlambda)
		gc()
		if(verbose) cat("done.\n")
	}
	fit$method = method
	rm(method,verbose)
	gc()
	fit$marker = "Successful"
	class(fit) = "scr"
	return(fit)
}

# Default printing function
print.scr = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.scr() has been terminated\n")
		return("Please refer to the manual")
	}
	if(x$method == "GACT") cat("This is a solution path using Graph Approximation via Correlation Thresholding (GACT) and length = ", length(x$path), "\n")
	if(!x$method == "GSS") cat("The dimension of prelected neighborhood after screening = ", nrow(x$ind.mat), "\n")
	cat("huge.scr() is an internal function, please refer to huge() and huge.select()")
}

# Default summary function
summary.scr = function(object, ...){
	if(object$marker == "Terminated"){
		cat("huge.scr() has been terminated\n")
		return("Please refer to the manual")
	}
	if(object$method == "GACT") cat("This is a solution path using Graph Approximation via Correlation Thresholding (GACT) and length = ", length(object$path), "\n")
	if(!object$method == "GSS") cat("The dimension of preselected neighborhood after screening = ", nrow(object$ind.mat), "\n")
	cat("huge.scr() is an internal function, please refer to huge() and huge.select()")
}

# Default plot function
plot.scr = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.scr() has been terminated\n")
		return("Please refer to the manual")
	}
	if(x$method == "GACT"){
		par(mfrow = c(1,1),pty = "s")
		plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "b",xlim = rev(range(x$lambda)))
	}
	if(!x$method == "GSS"){
		cat("No plot information avaiable for the Graph SURE Screening (GSS)","\n")
		cat("huge.scr() is an internal function. For more information, please refer to huge() and huge.select()")	}
}