#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.select(): (1)Rotation Information Criterion                      #
#                (2)Stability Approach to Regularization Selection      #
#                (3)Extended Bayesian Informaition Criterion            #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Apr 10th 2011                                                   #
# Version: 1.0.1                                                          #
#-----------------------------------------------------------------------#

## Main Function
huge.select = function(est, criterion = NULL, EBIC.gamma = 0.5, stars.thresh = 0.1, stars.subsample.ratio = NULL, stars.rep.num = 20, verbose = TRUE){

	gcinfo(FALSE)
		
	if(est$method == "MBGEL"&&is.null(criterion))
		criterion = "RIC"
	if(est$method == "GECT"&&is.null(criterion))
		criterion = "stars"
	if(est$method == "GLASSO"&&is.null(criterion))
		criterion = "RIC"
	
	n = nrow(est$data)
	d = ncol(est$data)
	nlambda = length(est$lambda)	
	
	if(criterion == "RIC")
	{
		if(verbose)
		{
			cat("Conducting Permutation Information Criterion (RIC) selection....")
        	flush.console()
        }
		
		out=.C("RIC",X = as.double(est$data),dd = as.integer(d),nn=as.integer(n),lambda_opt = as.double(0),PACKAGE="huge")
		est$opt.lambda = out$lambda_opt/n
		rm(out)
		gc()
		
		if(verbose){
			cat("done\n")
			flush.console()
		}
		
		if(verbose)
		{
			cat("Computing the optimal graph....")
			flush.console()
		}
		
		if(est$method == "MBGEL")
			est$refit = huge.MBGEL(est$data, lambda = est$opt.lambda, sym = est$sym, idx.mat = est$idx.mat, verbose = FALSE)$path[[1]]
		if(est$method == "GLASSO")
		{
			if(!is.null(est$w))
			{
				tmp = huge.glassoM(est$data, lambda = est$opt.lambda, cov.glasso = TRUE, verbose = FALSE)
				est$opt.w = tmp$w[[1]]
			}
			if(is.null(est$w))
				tmp = huge.glassoM(est$data, lambda = est$opt.lambda, verbose = FALSE)
			
			est$refit = tmp$path[[1]]
			est$opt.wi = tmp$wi[[1]]
			rm(tmp)
			gc()
		}
		if(est$method == "GECT")
			est$refit = huge.GECT(est$data, lambda = est$opt.lambda, verbose = FALSE)$path[[1]]
		
		est$opt.sparsity=sum(est$refit)/d/(d-1)
		
		if(verbose){
			cat("done\n")
			flush.console()
		}
	}
	
	if(criterion == "EBIC"&&est$method == "GLASSO")
	{
		est$EBIC.score = -2*est$loglik + log(n)*est$df + 2*EBIC.gamma*log(d*(d-1))*est$df
		est$opt.index = which.min(est$EBIC.score)
		est$refit = est$path[[est$opt.index]]
		est$opt.wi = est$wi[[est$opt.index]]
		if(!is.null(est$w))
			est$opt.w = est$w[[est$opt.index]]
  	est$opt.lambda = est$lambda[est$opt.index]
  	est$opt.sparsity = est$sparsity[est$opt.index]
	}		
	
	if(criterion == "stars"){
		if(is.null(stars.subsample.ratio))
		{
			if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
			if(n<=144) stars.subsample.ratio = 0.8
		} 
	
		est$merge = list()
		for(i in 1:nlambda) est$merge[[i]] = Matrix(0,d,d)
		
  		for(i in 1:stars.rep.num){
  			if(verbose){
				mes <- paste(c("Conducting Subsampling....in progress:", floor(100*i/stars.rep.num), "%"), collapse="")
   				cat(mes, "\r")
            	flush.console()	
			}
    		ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)
    		
    		if(est$method == "MBGEL")
    			tmp = huge.MBGEL(est$data[ind.sample,],lambda = est$lambda, scr = est$scr, idx.mat = est$idx.mat, sym = est$sym, verbose = FALSE)$path
       		if(est$method == "GECT")
    			tmp = huge.GECT(est$data[ind.sample,], lambda = est$lambda,verbose = FALSE)$path
    		if(est$method == "GLASSO")
    			tmp = huge.glassoM(est$data[ind.sample,], lambda = est$lambda, verbose = FALSE)$path
    			
    		for(i in 1:nlambda)
    			est$merge[[i]] = est$merge[[i]] + tmp[[i]]
    		
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
		for(i in 1:nlambda)
		{
			est$merge[[i]] = est$merge[[i]]/stars.rep.num
    		est$variability[i] = 4*sum(est$merge[[i]]*(1-est$merge[[i]]))/(d*(d-1))
    	}
       	est$opt.index = max(which.max(est$variability >= stars.thresh)[1]-1,1)
   		est$refit = est$path[[est$opt.index]]
  		est$opt.lambda = est$lambda[est$opt.index]
  		est$opt.sparsity = est$sparsity[est$opt.index]
  		if(est$method == "GLASSO")
  		{
  			est$opt.wi = est$wi[[est$opt.index]]
			  if(!is.null(est$w))
				  est$opt.w = est$w[[est$opt.index]]
  		}
	  }
    est$criterion = criterion
  	class(est) = "select"
    return(est)  	
}

#-----------------------------------------------------------------------#
# default printing function for class "select"                          #
#-----------------------------------------------------------------------#

print.select = function(x, ...)
{
	if(x$method == "GECT")
		cat("Model: Graph Estimation via Correlation Thresholding(GECT)\n")
	if(x$method == "GLASSO")
		cat("Model: Graphical Lasso (GLASSO)\n")
	if(x$method == "MBGEL")
		cat("Model: Meinshausen & Buhlmann Graph Estimation via Lasso(MBGEL)\n")

	cat("selection criterion:",x$criterion,"\n")
	if(x$NPN)
		cat("NonParaNormal transformed\n")
	if((x$method == "MBGEL")&&x$scr)
		cat("Graph SURE Screening (GSS): on\n")
	if(!is.null(x$theta))
		cat("True graph: included\n")
	if(is.null(x$theta))
		cat("True graph: not included\n")
	cat("Graph Dimension:",ncol(x$data),"\n")
	cat("sparsity level", x$opt.sparsity,"\n")
}

summary.select = function(object, ...){
	if(object$method == "GECT")
		cat("Model: Graph Estimation via Correlation Thresholding(GECT)\n")
	if(object$method == "GLASSO")
		cat("Model: Graphical Lasso (GLASSO)\n")
	if(object$method == "MBGEL")
		cat("Model: Meinshausen & Buhlmann Graph Estimation via Lasso(MBGEL)\n")

	cat("selection criterion:",object$criterion,"\n")
	if(object$NPN)
		cat("NonParaNormal transformed\n")
	if((object$method == "MBGEL")&&object$scr)
		cat("Graph SURE Screening (GSS): on\n")
	if(!is.null(object$theta))
		cat("True graph: included\n")
	if(is.null(object$theta))
		cat("True graph: not included\n")
	cat("Graph Dimension:",ncol(object$data),"\n")
	cat("sparsity level", object$opt.sparsity,"\n")
}

plot.select = function(x, ...){
	par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
    
    g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
	layout.grid = layout.fruchterman.reingold(g)
	
	plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA)	  
    plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Solution path sparsity levels")
    lines(x$opt.lambda,x$opt.sparsity,type = "p")  
}