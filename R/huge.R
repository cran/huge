#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge(): High-dimensional Undirected Graph Estimation via              #
#		  (1) Lasso 	   		                                        #
#		  (2) Correlation Thresholding                                  #
#		  (3) Graphical Lasso            	                            #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Feb 28th 2011                                                   #
# Version: 1.0                                                          #
#-----------------------------------------------------------------------#

## Main function
huge = function(L, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, NPN = FALSE, NPN.func = "shrinkage", NPN.thresh = NULL, method = "MBGEL", scr = NULL, scr.num = NULL, cov.glasso = FALSE, sym = "or", verbose = TRUE)
{	
	gcinfo(FALSE)
	est = list()
	est$method = method
	
	if(is.list(L))
	{
		n = nrow(L$data)
		d = ncol(L$data)
		est$data = L$data
		if(!is.null(L$theta))	est$theta = L$theta
	}
	
	if(!is.list(L))
	{
		n = nrow(L)
		d = ncol(L)
		est$data = L
	}
	rm(L)
	gc()		
	
	# Nonparanormal transformation
	if(NPN)
	{
		est$data = huge.NPN(est$data,NPN.func = NPN.func, NPN.thresh = NPN.thresh, verbose = verbose)$data
		rm(NPN.thresh,NPN.func)
		gc()
	}
	est$NPN = NPN
	
	if(method == "GECT")
	{
		fit = huge.GECT(est$data, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		rm(fit)
		gc()
	}
	
	if(method == "MBGEL")
	{	
		fit = huge.MBGEL(est$data, lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, scr = scr, scr.num = scr.num, sym = sym, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		est$rss = fit$rss
		est$df = fit$df
		est$idx_mat = fit$idx_mat
		est$sym = sym
		est$scr = fit$scr
		rm(fit,sym)
		gc()
	}
	
	
	if(method == "GLASSO")
	{
		fit = huge.glassoM(est$data, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, cov.glasso = cov.glasso, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$wi = fit$wi
		est$df = fit$df
		est$sparsity = fit$sparsity
		est$loglik = fit$loglik
		if(cov.glasso)
			est$w = fit$w
		rm(fit)
		gc()
	}			
			
	rm(n,d,NPN,scr,lambda,lambda.min.ratio,nlambda,verbose)
	gc()
	
	class(est) = "huge"
	return(est)
}

print.huge = function(x, ...)
{	
	if(x$method == "GECT")
		cat("Model: Graph Approximation via Correlation Thresholding (GECT)\n")
	if(x$method == "GLASSO")
		cat("Model: Graphical Lasso (GLASSO)\n")
	if(x$method == "MBGEL")
		cat("Model: Meinshausen & Buhlmann Graph Estimation via Lasso (MBGEL)\n")
	
	if(x$NPN) cat("Nonparanormal (NPN) transformation: on\n")
	
	if((x$method == "MBGEL")&&(x$scr)) cat("Graph SURE Screening (GSS): on\n")

	if(is.null(x$theta)) cat("True graph: not included\n")
	if(!is.null(x$theta)) cat("True graph: included\n")
	
	cat("Path length:",length(x$lambda),"\n")
	cat("Graph Dimension:",ncol(x$data),"\n")
	cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
}

summary.huge = function(object, ...)
{	
	if(object$method == "GECT")
		cat("Model: Graph Estimation via Correlation Thresholding (GECT)\n")
	if(object$method == "GLASSO")
		cat("Model: Graphical Lasso (GLASSO)\n")
	if(object$method == "MBGEL")
		cat("Model: Meinshausen & Buhlmann Graph Estimation via Lasso (MBGEL)\n")
	
	if(object$NPN) cat("Nonparanormal (NPN) transformation: on\n")
	
	if((object$method == "MBGEL")&&object$scr) cat("Graph SURE Screening (GSS): on\n")

	if(is.null(object$theta)) cat("True graph: not included\n")
	if(!is.null(object$theta)) cat("True graph: included\n")
	
	cat("Path length:",length(object$lambda),"\n")
	cat("Graph Dimension:",ncol(object$data),"\n")
	cat("Sparsity level:",min(object$sparsity),"----->",max(object$sparsity),"\n")

}

plot.huge = function(x, align = FALSE, ...){
	gcinfo(FALSE)
	
	if(length(x$lambda) == 1)	par(mfrow = c(1, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	if(length(x$lambda) == 2)	par(mfrow = c(1, 3), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	if(length(x$lambda) >= 3)	par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	
	if(length(x$lambda) <= 3)	z.final = 1:length(x$lambda)
	
	if(length(x$lambda) >=4){
		z.max = max(x$sparsity)
		z.min = min(x$sparsity)
		z = z.max - z.min
		z.unique = unique(c(which(x$sparsity>=(z.min + 0.03*z))[1],which(x$sparsity>=(z.min + 0.07*z))[1],which(x$sparsity>=(z.min + 0.15*z))[1]))

		
		if(length(z.unique) == 1){
			if(z.unique<(length(x$lambda)-1))	z.final = c(z.unique,z.unique+1,z.unique+2)
			if(z.unique==(length(x$lambda)-1)) z.final = c(z.unique-1,z.unique,z.unique+1)
			if(z.unique==length(x$lambda)) 	z.final = c(z.unique-2,z.unique-1,z.unique)
		}
		
		if(length(z.unique) == 2){
			if(diff(z.unique)==1){
				if(z.unique[2]<length(x$lambda)) z.final = c(z.unique,z.unique[2]+1) 
				if(z.unique[2]==length(x$lambda)) z.final = c(z.unique[1]-1,z.unique)
			}
			if(diff(z.unique)>1) z.final = c(z.unique[1],z.unique[1]+1,z.unique[2])
		}
		
		if(length(z.unique) == 3) z.final = z.unique
		
		rm(z.max,z.min,z,z.unique)
		gc()
		
	}
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Sparsity vs. Regularization")
	
	lines(x$lambda[z.final],x$sparsity[z.final],type = "p")
	
	if(align){
		layout.grid = layout.fruchterman.reingold(graph.adjacency(as.matrix(x$path[[z.final[length(z.final)]]]), mode="undirected", diag=FALSE))
		for(i in z.final){
			g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g)
			gc()
		}
	rm(layout.grid)
	}
	if(!align){
		for(i in z.final){
			g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
			layout.grid = layout.fruchterman.reingold(g)
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g,layout.grid)
			gc()
		}
	}
	if(align) cat("Three plotted graphs are aligned according to the third graph\n")
}