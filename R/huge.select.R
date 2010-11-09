#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.select(): StARS selection and Extended BIC                       #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@gmail.com>; <hanliu@cs.jhu.edu>                     #
# Date: Nov 9th, 2010                                                   #
# Version: 0.7                                                          #
#-----------------------------------------------------------------------#

huge.stars = function(R.path, stars.thresh = 0.1, verbose = TRUE){
		
	if(verbose) cat('Statbility selection....\n')
	n.subs = length(R.path)
  	n.lambda = length(R.path[[1]])
  	d = ncol(R.path[[1]][[1]])

	D = rep(0,n.lambda)
	for(i in 1:n.lambda){
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
			if(i==1) cat("The threshold is too small! Try a different set of regularization/thresholding parameters....\n")
			rm(n.subs,n.lambda,d, R.path,D)
   			gc(gcinfo(verbose = FALSE))
   			opt.index = max(c(1,i-1))
			return(opt.index)
		}
	}
	rm(n.subs,n.lambda,d,R.path)
   	gc(gcinfo(verbose = FALSE))
   	return(i)
}

#-----------------------------------------------------------------------#
# model selection function for class "huge"                             #
#-----------------------------------------------------------------------#
huge.select = function(est, criterion = NULL, EBIC.gamma = 0.5, stars.thresh = 0.1, sample.ratio = NULL, rep.num = 20, verbose = TRUE){

	if(est$marker == "Terminated"){
		cat("huge.select() will be terminated....\n")
		class(est) = "select"
		return(est)
	}
		
	if(!est$approx&&is.null(criterion)) criterion = "EBIC"
	if(est$approx&&is.null(criterion)) criterion = "stars"
	n = nrow(est$data)
	d = ncol(est$data)
	k = length(est$ind.group)
	n.lambda = length(est$lambda)
	
	if(est$type == "subgraph solution path") est$type = "optimal subgraph"
	if(est$type == "fullgraph solution path") est$type = "optimal fullgraph"	
	if(criterion == "EBIC"){
		scr.num = min(floor(2*n/3),d-1)
		ind.mat = huge.scr(est$data, ind.group = est$ind.group, scr.num  = scr.num, verbose = FALSE)$ind.mat
		EBIC.ind = matrix(0,k,n.lambda)
		for(i in est$ind.group){
			sigma2hat = sum((est$data[,est$ind.group[i]] - est$data[,ind.mat[,i]]%*%lm.ridge(est$data[,est$ind.group[i]]~est$data[,ind.mat[,i]])$coef)^2)/(n-scr.num-1)
			EBIC.ind[i,] = est$rss[i,]/n/sigma2hat + log(n)*est$df[i,]/n + 2*est$df[i,]*EBIC.gamma*log(d)/n
		}	
		
		est$EBIC.score = apply(EBIC.ind,2,sum)
		est$opt.index = which.min(est$EBIC.score)
	}
	
		
	if(criterion == "stars"){
		if(is.null(sample.ratio)){
			if(est$n>144) sample.ratio = 10*sqrt(n)/n
			if(est$n<=144) sample.ratio = 0.8
		} 
	
		R.path = list()
  		for(i in 1:rep.num){
  			if(verbose) cat('subsampling for stars...',i,'\n')
    		ind.sample = sample(c(1:n), floor(n*sample.ratio), replace=FALSE)
    		
    		if(!est$approx){
    			if(is.null(est$ind.mat))
    				R.path[[i]] = huge.subgraph(est$data[ind.sample,],ind.group = est$ind.group, alpha = est$alpha, lambda = est$lambda, sym = est$sym, verbose = verbose)$path
    			if(!is.null(est$ind.mat))
    				R.path[[i]] = huge.subgraph(est$data[ind.sample,],ind.group = est$ind.group, ind.mat = est$ind.mat, alpha = est$alpha, lambda = est$lambda, sym = est$sym, verbose = verbose)$path
    		}
    		
    		if(est$approx)
    			R.path[[i]] = huge.scr(est$data[ind.sample,],ind.group = est$ind.group,lambda = est$lambda, approx = TRUE)$path

    		rm(ind.sample)
   			gc(gcinfo(verbose = FALSE))
		}

  		est$opt.index = huge.stars(R.path, stars.thresh, verbose)
  		est$merge = Matrix(0,k,k)
  		for(i in 1:rep.num)
  			est$merge = est$merge + R.path[[i]][[est$opt.index]]
  		rm(R.path)
   		gc(gcinfo(verbose = FALSE))
   	}
  	
  	est$refit = est$path[[est$opt.index]]
  	est$opt.lambda = est$lambda[est$opt.index]
  	est$opt.sparsity = est$sparsity[est$opt.index]
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
	if(x$approx) cat("Model: correlation graph estimation\n")
	
	if(!x$approx){
		if(x$alpha <1) cat("Model: Meinshausen &Buhlmann graph estimation using elastic net --->",x$type,"\n")
		if(x$alpha ==1) cat("Model: Meinshausen &Buhlmann graph estimation using lasso --->",x$type,"\n")
	}
	cat("selection criterion:",x$criterion,"\n")
	if(x$npn) cat("nonparanormal transformed\n")
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
	if(object$approx) cat("Model: correlation graph estimation\n")
	
	if(!object$approx){
		if(object$alpha <1) cat("Model: Meinshausen &Buhlmann graph estimation using elastic net --->",object$type,"\n")
		if(object$alpha ==1) cat("Model: Meinshausen &Buhlmann graph estimation using lasso --->",object$type,"\n")
	}
	cat("selection criterion:",object$criterion,"\n")
	if(object$npn) cat("nonparanormal transformed\n")
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