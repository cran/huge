#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge(): Onestep estimation for high-dimensional undirected graph      #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 21st 2010                                                   #
# Version: 0.9                                                          #
#-----------------------------------------------------------------------#

## Main function
huge = function(L, ind.group = NULL, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.1, alpha = 1, sym = "or", npn = TRUE, npn.func = "shrinkage", npn.thresh = NULL, method = "GEL", scr = TRUE, scr.num = NULL, verbose = TRUE){	
	gcinfo(FALSE)
	est = list()
	est$method = method
	
	if(is.list(L)){
		n = nrow(L$data)
		d = ncol(L$data)
		if(d < 3){
			cat("The fullgraph dimension < 3 and huge() will be teminated....\n")
			cat("Please refer to Pearson's product-moment correlation....\n")
			rm(L,alpha,sym,npn.func,method,scr,verbose)
			gc()
			est$marker = "Terminated"
			class(est) = "huge"
			return(est)
		}
		if(!is.null(scr.num))
			if(scr.num == 1){
				cat("The neighborhood size < 2 and huge() will be teminated.\n")
				cat("Please refer to Pearson's product-moment correlation....\n")
				rm(L,alpha,sym,npn.func,method,scr,verbose)
				gc()
				est$marker = "Terminated"
				class(est) = "huge"
				return(est)
		}
		est$data = L$data
		if(!is.null(L$theta))	est$theta = L$theta
	}
	
	if(!is.list(L)){
		n = nrow(L)
		d = ncol(L)
		if(d < 3){
			cat("The fullgraph dimension < 3 and huge() will be teminated....\n")
			cat("Please refer to Pearson's product-moment correlation....\n")
			rm(L,alpha,sym,npn.func,method,scr,verbose)
			gc()
			est$marker = "Terminated"
			class(est) = "huge"
			return(est)
		}
		if(!is.null(scr.num))
			if(scr.num == 1){
				cat("The neighborhood size < 2 and huge() will be teminated.\n")
				cat("Please refer to Pearson's product-moment correlation....\n")
				est$marker = "Terminated"
				rm(L,alpha,sym,npn.func,method,scr,verbose)
				gc()
				class(est) = "huge"
				return(est)
		}
		est$data = L
	}
	rm(L)
	gc()	
	
	if(is.null(ind.group))	ind.group = c(1:d)
	k = length(ind.group)	
	
	# Nonparanormal transformation
	if(npn){
		est$data = huge.npn(est$data,npn.func = npn.func, npn.thresh = npn.thresh, verbose = verbose)$data
		rm(npn.thresh,npn.func)
		gc()
	}
	
	if(method == "GLASSO"){
		if(is.null(lambda)){
			if(is.null(nlambda)) nlambda = 10
		}
		fit = huge.glassoM(est$data, ind.group = ind.group, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$wi = fit$wi
		est$df = fit$df
		est$sparsity = fit$sparsity
		est$loglik = fit$loglik
		rm(fit)
		gc()
	}			
			
	if(method == "GACT"){
		if(is.null(lambda)){
			if(is.null(nlambda)) nlambda = 30
		}
		fit = huge.scr(est$data, ind.group = ind.group, method = method, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		rm(fit)
		gc()
	}
	
	if(method == "GEL"){		
		if(!scr){
			ind.mat = NULL
			est$scr = FALSE
		}
		if(scr){
			if(is.null(scr.num)&&(n >= d)){
				ind.mat = NULL
				est$scr = FALSE
			}
			
			if(is.null(scr.num)&&(n < d)){
				scr.num = n-1
				est$scr = FALSE
			}
			
			if(!is.null(scr.num)){
				if(scr.num >= (d-1)){
					if(verbose) cat("The specified scr.num >= (d-1) and Graph SURE Screening will be skipped.\n")
					ind.mat = NULL
					est$scr = FALSE
				}
				if(scr.num < (d-1)){
					ind.mat = huge.scr(est$data, ind.group = ind.group, scr.num = scr.num, verbose = verbose)$ind.mat
					est$scr = TRUE
				}
			}
		}
		if(is.null(lambda)){
			if(is.null(nlambda)) nlambda = 10
		}
		fit = huge.subgraph(est$data, ind.group = ind.group, ind.mat = ind.mat, alpha = alpha, lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, sym = sym, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		est$rss = fit$rss
		est$df = fit$df
		rm(fit)
		gc()
		
		est$ind.mat = ind.mat
		est$alpha = alpha
		est$sym = sym
		rm(ind.mat,alpha,sym,method)
		gc()
	}
	est$ind.group = ind.group
	if(k<d) est$type = "Subgraph solution path"
	if(k==d) est$type = "Fullgraph solution path"
	est$npn = npn
	
	rm(ind.group,n,k,d,npn,scr,lambda,lambda.min.ratio,nlambda,verbose)
	gc()	
	
	est$marker = "Successful"
	
	class(est) = "huge"
	return(est)
}

print.huge = function(x, ...){
	
	if(x$marker == "Terminated"){
		cat("huge() has been terminated\n")
		return("Please refer to the manual")
	}
	
	if(x$method == "GACT") cat("Model: Graph Approximation via Correlation Thresholding (GACT) --->",x$type,"\n")
	if(x$method == "GLASSO") cat("Model: Graphical Lasso (GLASSO) --->",x$type,"\n")
	if(x$method == "GEL"){
		if(x$alpha <1) cat("Model: Meinshausen &Buhlmann Graph Estimation via Elastic Net --->",x$type,"\n")
		if(x$alpha ==1) cat("Model: Meinshausen &Buhlmann Graph Estimation via Lasso (GEL) --->",x$type,"\n")
	}
	
	if(x$npn) cat("Nonparanormal (NPN) transformation: on\n")
	
	if((x$method == "GEL")&&x$scr) cat("Graph SURE Screening (GSS): on\n")

	if(is.null(x$theta)) cat("True graph: not included\n")
	if(!is.null(x$theta)) cat("True graph: included\n")
	
	cat("Path length:",length(x$lambda),"\n")
	cat("Subgraph dimension:",length(x$ind.group),"\n")
	cat("Fullgraph dimension:",ncol(x$data),"\n")
	cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
}

summary.huge = function(object, ...){
	
	if(object$marker == "Terminated"){
		cat("huge() has been terminated\n")
		return("Please refer to the manual")
	}
	if(object$method == "GACT") cat("Model: Graph Approximation via Correlation Thresholding (GACT) --->",object$type,"\n")
	if(object$method == "GLASSO") cat("Model: Graphical Lasso (GLASSO) --->",object$type,"\n")
	if(object$method == "GEL"){
		if(object$alpha <1) cat("Model: Meinshausen &Buhlmann Graph Estimation via Elastic Net --->",object$type,"\n")
		if(object$alpha ==1) cat("Model: Meinshausen &Buhlmann Graph Estimation via Lasso (GEL)--->",object$type,"\n")
	}
	
	if(object$npn) cat("Nonparanormal (NPN) transformation: on\n")
	if((object$method == "GEL")&&object$scr) cat("Graph SURE Screening (GSS): on\n")

	if(is.null(object$theta)) cat("True graph: not included\n")
	if(!is.null(object$theta)) cat("True graph: included\n")
	
	cat("Path length:",length(object$lambda),"\n")
	cat("Subgraph dimension:",length(object$ind.group),"\n")
	cat("Fullgraph dimension:",ncol(object$data),"\n")
	cat("Sparsity level:",min(object$sparsity),"----->",max(object$sparsity),"\n")
}

plot.huge = function(x, align = FALSE, ...){
	gcinfo(FALSE)
	if(x$marker == "Terminated"){
		cat("huge() has been terminated\n")
		return("Please refer to the manual")
	}
	
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
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g)
			gc()
		}
	rm(layout.grid)
	}
	if(!align){
		for(i in z.final){
			g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
			layout.grid = layout.fruchterman.reingold(g)
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g,layout.grid)
			gc()
		}
	}
	if(align) cat("Three plotted graphs are aligned according to the third graph\n")
}