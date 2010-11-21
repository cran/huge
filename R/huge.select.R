#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.select(): (1)Permutation Information Criterion                   #
#                (2)Stability Approach to Regularization Selection      #
#                (3)Extended Bayesian Informaition Criterion            #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 21st 2010                                                   #
# Version: 0.9                                                          #
#-----------------------------------------------------------------------#

## Main Function
huge.select = function(est, criterion = NULL, r.num = 200, EBIC.gamma = 0.5, stars.thresh = 0.1, stars.subsample.ratio = NULL, stars.rep.num = 20, verbose = TRUE){

	gcinfo(FALSE)
	
	if(est$marker == "Terminated"){
		cat("huge.select() will be terminated....\n")
		class(est) = "select"
		return(est)
	}
		
	if(est$method == "GEL"&&is.null(criterion)) criterion = "PIC"
	if(est$method == "GACT"&&is.null(criterion)) criterion = "stars"
	if(est$method == "GLASSO"&&is.null(criterion)) criterion = "EBIC"
	
	n = nrow(est$data)
	d = ncol(est$data)
	ind.group = est$ind.group
	k = length(est$ind.group)
	nlambda = length(est$lambda)
	
	if(est$type == "Subgraph solution path") est$type = "Optimal subgraph"
	if(est$type == "Fullgraph solution path") est$type = "Optimal fullgraph"	
	
	if(criterion == "PIC"){
		lambda.tmp = rep(0,k)
		if(est$scr)
			for(i in 1:k)
				lambda.tmp[i] = max(abs(t(est$data[,ind.group[i]]%*%est$data[,est$ind.mat[,ind.group[i]]])))
		if(!est$scr)
			for(i in 1:k)
				lambda.tmp[i] = max(abs(t(est$data[,ind.group[i]]%*%est$data[,-ind.group[i]])))
		lambda.max = max(lambda.tmp)
		
		lambda.all = rep(0,r.num)
		for(r in 1:r.num){
			if(verbose){
				mes <- paste(c("Conducting Permutation Information Criterion (PIC) selection....in progress:", floor(100*r/r.num), "%"), collapse="")
   				cat(mes, "\r")
            	flush.console()	
			}
			if(est$scr){
				for(i in 1:k){
					rperm = sample(1:n,n)
					lambda.tmp[i] = max(abs(t(est$data[rperm,ind.group[i]]%*%est$data[,est$ind.mat[,ind.group[i]]])))
				}
				lambda.all[r] = max(lambda.tmp)
			}
			if(!est$scr){
				for(i in 1:k){
					rperm = sample(1:n,n)
					lambda.tmp[i] = max(abs(t(est$data[rperm,ind.group[i]]%*%est$data[,-ind.group[i]])))
				}
			lambda.all[r] = max(lambda.tmp)
			}
		}
		if(verbose){
			mes = "Conducting Permutation Information Criterion (PIC) selection....done.                 "
        	cat(mes, "\r")
        	cat("\n")
        	flush.console()
        }
		lambda.min = min(lambda.all)
		lambda = exp(seq(log(lambda.max),log(lambda.min),length = 5))/n/est$alpha
		est$opt.lambda = lambda.min/n/est$alpha
		if(verbose) cat("Computing the optimal graph\n")
		est$refit = huge.subgraph(est$data, ind.group = est$ind.group, ind.mat = est$ind.mat, alpha = est$alpha,lambda = lambda,sym = est$sym, verbose = verbose)$path[[5]]
		est$opt.sparsity=sum(est$refit)/k/(k-1)
	}
	
	if(criterion == "EBIC"&&est$method == "GEL"){
		scr.num = min(floor(2*n/3),d-1)
		ind.mat = huge.scr(est$data, ind.group = est$ind.group, scr.num  = scr.num, verbose = FALSE)$ind.mat
		EBIC.ind = matrix(0,k,nlambda)
		for(i in 1:k){
			if(verbose){
				mes <- paste(c("Conducting Extended Bayesian Information Criterion (EBIC) selection....in progress:", floor(100*i/k), "%"), collapse="")
   				cat(mes, "\r")
            	flush.console()	
			}
			sigma2hat = sum((est$data[,est$ind.group[i]] - est$data[,ind.mat[,est$ind.group[i]]]%*%lm.ridge(est$data[,est$ind.group[i]]~est$data[,ind.mat[,est$ind.group[i]]])$coef)^2)/(n-scr.num-1)
			EBIC.ind[i,] = est$rss[i,]/n/sigma2hat + log(n)*est$df[i,]/n + 2*est$df[i,]*EBIC.gamma*log(d)/n
		}
		if(verbose){
			mes = "Conducting Extended Bayesian Information Criterion (EBIC) selection....done.                 "
        	cat(mes, "\r")
        	cat("\n")
        	flush.console()
        }	
		
		est$EBIC.score = apply(EBIC.ind,2,sum)
		rm(EBIC.ind)
		gc()
		est$opt.index = which.min(est$EBIC.score)
		est$refit = est$path[[est$opt.index]]
  		est$opt.lambda = est$lambda[est$opt.index]
  		est$opt.sparsity = est$sparsity[est$opt.index]
	}
	
	if(criterion == "EBIC"&&est$method == "GLASSO"){
		est$EBIC.score = -2*est$loglik + log(n)*est$df + 2*EBIC.gamma*log(d*(d-1))*est$df
		est$opt.index = which.min(est$EBIC.score)
		est$refit = est$path[[est$opt.index]]
		est$opt.wi = est$wi[[est$opt.index]]
  		est$opt.lambda = est$lambda[est$opt.index]
  		est$opt.sparsity = est$sparsity[est$opt.index]
	}		
	if(criterion == "stars"){
		if(is.null(stars.subsample.ratio)){
			if(est$n>144) stars.subsample.ratio = 10*sqrt(n)/n
			if(est$n<=144) stars.subsample.ratio = 0.8
		} 
	
		est$merge = list()
		for(i in 1:nlambda) est$merge[[i]] = Matrix(0,k,k)
		
  		for(i in 1:stars.rep.num){
  			if(verbose){
				mes <- paste(c("Conducting Subsampling....in progress:", floor(100*i/stars.rep.num), "%"), collapse="")
   				cat(mes, "\r")
            	flush.console()	
			}
    		ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)
    		
    		if(est$method == "GEL"){
    			if(is.null(est$ind.mat))
    				tmp = huge.subgraph(est$data[ind.sample,],ind.group = est$ind.group, alpha = est$alpha, lambda = est$lambda, sym = est$sym, verbose = FALSE)$path
    			if(!is.null(est$ind.mat))
    				tmp = huge.subgraph(est$data[ind.sample,],ind.group = est$ind.group, ind.mat = est$ind.mat, alpha = est$alpha, lambda = est$lambda, sym = est$sym, verbose = FALSE)$path
    		}
       		if(est$method == "GACT")
    			tmp = huge.scr(est$data[ind.sample,],ind.group = est$ind.group,lambda = est$lambda, method = "GACT",verbose = FALSE)$path
    		if(est$method == "GLASSO")
    			tmp = huge.glassoM(est$data[ind.sample,],ind.group = est$ind.group, lambda = est$lambda, verbose = FALSE)$path
    			
    		for(i in 1:nlambda)	est$merge[[i]] = est$merge[[i]] + tmp[[i]]
    		rm(ind.sample,tmp)
   			gc()
		}
		
		if(verbose){
			mes = "Conducting Subsampling....done.                 "
        	cat(mes, "\r")
        	cat("\n")
        	flush.console()
        }
        
        est$variability = rep(0,nlambda)
		for(i in 1:nlambda){
			est$merge[[i]] = est$merge[[i]]/stars.rep.num
    		est$variability[i] = 4*sum(est$merge[[i]]*(1-est$merge[[i]]))/(k*(k-1))
    	}
       	est$opt.index = max(which.max(est$variability >= stars.thresh)[1]-1,1)
   		est$refit = est$path[[est$opt.index]]
  		est$opt.lambda = est$lambda[est$opt.index]
  		est$opt.sparsity = est$sparsity[est$opt.index]
	}
    est$criterion = criterion
  	class(est) = "select"
    return(est)  	
}

#-----------------------------------------------------------------------#
# default printing function for class "select"                          #
#-----------------------------------------------------------------------#

print.select = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.select() has been terminated\n")
		return("Please refer to the manual")
	}
	if(x$method == "GACT") cat("Model: Graph Approximation via Correlation Thresholding(GACT)\n")
	if(x$method == "GLASSO") cat("Model: Graphical Lasso (GLASSO)\n")
	if(x$method == "GEL"){
		if(x$alpha <1) cat("Model: Meinshausen & Buhlmann Graph Estimation using Elastic Net --->",x$type,"\n")
		if(x$alpha ==1) cat("Model: Meinshausen & Buhlmann Graph Estimation using Lasso(GEL) --->",x$type,"\n")
	}
	cat("selection criterion:",x$criterion,"\n")
	if(x$npn) cat("nonparanormal transformed\n")
	if(x$method == "GEL")
		if(x$scr) cat("graph screened\n")
	if(!is.null(x$theta)) cat("ground truth theta: included\n")
	if(is.null(x$theta)) cat("ground truth theta: not included\n")
	cat("subgraph dimension:",length(x$ind.group),"\n")
	cat("fullgraph dimension:",ncol(x$data),"\n")
	cat("sparsity level", x$opt.sparsity,"\n")
}

summary.select = function(object, ...){
	if(object$marker == "Terminated"){
		cat("huge.select() has been terminated\n")
		return("Please refer to the manual")
	}
	if(object$method == "GACT") cat("Model: Graph Approximation via Correlation Thresholding(GACT)\n")
	if(object$method == "GLASSO") cat("Model: Graphical Lasso (GLASSO)\n")
	if(object$method == "GEL"){
		if(object$alpha <1) cat("Model: Meinshausen & Buhlmann Graph Estimation using Elastic Net --->",object$type,"\n")
		if(object$alpha ==1) cat("Model: Meinshausen & Buhlmann Graph Estimation using Lasso (GEL) --->",object$type,"\n")
	}
	cat("selection criterion:",object$criterion,"\n")
	if(object$npn) cat("nonparanormal transformed\n")
	if(object$method == "GEL")
		if(object$scr) cat("graph screened\n")
	if(!is.null(object$theta)) cat("ground truth theta: included\n")
	if(is.null(object$theta)) cat("ground truth theta: not included\n")
	cat("subgraph dimension:",length(object$ind.group),"\n")
	cat("fullgraph dimension:",ncol(object$data),"\n")
	cat("sparsity level", object$opt.sparsity,"\n")
}

plot.select = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.select() has been terminated\n")
		return("Please refer to the manual")
	}
	par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
    
    g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
	layout.grid = layout.fruchterman.reingold(g)
	
	plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA)	  
    plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Solution path sparsity levels")
    lines(x$opt.lambda,x$opt.sparsity,type = "p")  
}