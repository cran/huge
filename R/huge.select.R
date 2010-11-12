#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.select(): (1)Permutation Information Criterion                   #
#                (2)Stability Approach to Regularization Selection      #
#                (3)Extended Bayesian Informaition Criterion            #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 12th, 2010                                                  #
# Version: 0.8                                                          #
#-----------------------------------------------------------------------#

## Main Function
huge.select = function(est, criterion = NULL, r.num = 200, EBIC.gamma = 0.5, stars.thresh = 0.1, stars.subsample.ratio = NULL, stars.rep.num = 20, verbose = TRUE){

	if(est$marker == "Terminated"){
		cat("huge.select() will be terminated....\n")
		class(est) = "select"
		return(est)
	}
		
	if(!est$approx&&is.null(criterion)) criterion = "PIC"
	if(est$approx&&is.null(criterion)) criterion = "stars"
	
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
		if(verbose) cat("\n")
		if(verbose) cat("Computing the optimal graph\n")
		est$refit = huge.subgraph(est$data,est$ind.group,est$ind.mat,est$alpha,lambda,sym = est$sym,verbose = verbose)$path[[5]]
		est$opt.sparsity=sum(est$refit)/k/(k-1)
	}
	
	if(criterion == "EBIC"){
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
		est$opt.index = which.min(est$EBIC.score)
		est$refit = est$path[[est$opt.index]]
  		est$opt.lambda = est$lambda[est$opt.index]
  		est$opt.sparsity = est$sparsity[est$opt.index]
	}
	
		
	if(criterion == "stars"){
		if(is.null(stars.subsample.ratio)){
			if(est$n>144) stars.subsample.ratio = 10*sqrt(n)/n
			if(est$n<=144) stars.subsample.ratio = 0.8
		} 
	
		R.path = list()
  		for(i in 1:stars.rep.num){
  			if(verbose){
				mes <- paste(c("Conducting Subsampling....in progress:", floor(100*i/stars.rep.num), "%"), collapse="")
   				cat(mes, "\r")
            	flush.console()	
			}
    		ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)
    		
    		if(!est$approx){
    			if(is.null(est$ind.mat))
    				R.path[[i]] = huge.subgraph(est$data[ind.sample,],ind.group = est$ind.group, alpha = est$alpha, lambda = est$lambda, sym = est$sym, verbose = FALSE)$path
    			if(!is.null(est$ind.mat))
    				R.path[[i]] = huge.subgraph(est$data[ind.sample,],ind.group = est$ind.group, ind.mat = est$ind.mat, alpha = est$alpha, lambda = est$lambda, sym = est$sym, verbose = FALSE)$path
    		}
       		if(est$approx)
    			R.path[[i]] = huge.scr(est$data[ind.sample,],ind.group = est$ind.group,lambda = est$lambda, approx = TRUE,verbose = FALSE)$path

    		rm(ind.sample)
   			gc(gcinfo(verbose = FALSE))
		}
		if(verbose){
			mes = "Conducting Subsampling....done.                 "
        	cat(mes, "\r")
        	cat("\n")
        	flush.console()
        }

  		est$opt.index = huge.stars(R.path, stars.thresh, verbose)
  		est$merge = Matrix(0,k,k)
  		for(i in 1:stars.rep.num)
  			est$merge = est$merge + R.path[[i]][[est$opt.index]]
  		rm(R.path)
   		gc(gcinfo(verbose = FALSE))
   		est$refit = est$path[[est$opt.index]]
  		est$opt.lambda = est$lambda[est$opt.index]
  		est$opt.sparsity = est$sparsity[est$opt.index]
	}
    est$criterion = criterion
  	class(est) = "select"
    return(est)  	
}




huge.stars = function(R.path, stars.thresh = 0.1, verbose = TRUE){
		
	if(verbose) cat('Conducting Statbility Approach to Regularization Selection (StARS)....')
	n.subs = length(R.path)
  	nlambda = length(R.path[[1]])
  	d = ncol(R.path[[1]][[1]])

	D = rep(0,nlambda)
	for(i in 1:nlambda){
		P = matrix(0, d, d)
		for(j in 1:n.subs){
			P = P + R.path[[j]][[i]]
		}
		P = P/n.subs
		D[i] = 4*sum(P*(1-P))/(d*(d-1))
		rm(P)
   		gc(gcinfo(verbose = FALSE))
		D[i] = max(D[1:i])
		if(D[i]>stars.thresh){
			if(i==1) cat("\n The threshold is too small! Try a different set of regularization/thresholding parameters....\n")
			rm(n.subs,nlambda,d, R.path,D)
   			gc(gcinfo(verbose = FALSE))
   			opt.index = max(c(1,i-1))
   			if(verbose) cat("done.\n")
			return(opt.index)
		}
	}
	rm(n.subs,nlambda,d,R.path)
   	gc(gcinfo(verbose = FALSE))
   	if(verbose) cat("done.\n")
   	return(i)
}

#-----------------------------------------------------------------------#
# default printing function for class "select"                          #
#-----------------------------------------------------------------------#

print.select = function(x, ...){
	if(x$marker == "Terminated"){
		cat("huge.select() has been terminated\n")
		return("Please refer to the manual")
	}
	if(x$approx) cat("Model: Graph Approximation via Correlation Thresholding(GACT)\n")
	
	if(!x$approx){
		if(x$alpha <1) cat("Model: Meinshausen & Buhlmann Graph Estimation using Elastic Net --->",x$type,"\n")
		if(x$alpha ==1) cat("Model: Meinshausen & Buhlmann Graph Estimation using Lasso(GEL) --->",x$type,"\n")
	}
	cat("selection criterion:",x$criterion,"\n")
	if(x$npn) cat("nonparanormal transformed\n")
	if(!x$approx)
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
	if(object$approx) cat("Model: Graph Approximation via Correlation Thresholding (GACT)\n")

	
	if(!object$approx){
		if(object$alpha <1) cat("Model: Meinshausen & Buhlmann Graph Estimation using Elastic Net --->",object$type,"\n")
		if(object$alpha ==1) cat("Model: Meinshausen & Buhlmann Graph Estimation using Lasso (GEL) --->",object$type,"\n")
	}
	cat("selection criterion:",object$criterion,"\n")
	if(object$npn) cat("nonparanormal transformed\n")
	if(!object$approx)
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
	if(sum(x$refit) == 0) cat("The optimal graph is a null graph.")
	if(sum(x$refit) > 0){
	par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
    
    g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
	layout.grid = layout.fruchterman.reingold(g)
	
	plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA)	  
    plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Solution path sparsity levels")
    lines(x$opt.lambda,x$opt.sparsity,type = "p")  
    }
}