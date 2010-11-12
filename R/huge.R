#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge(): Onestep estimation for high-dimensional undirected graph      #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 12th 2010                                                   #
# Version: 0.8                                                          #
#-----------------------------------------------------------------------#

## Main function
huge = function(L, ind.group = NULL, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, alpha = 1, sym = "or", npn = TRUE, npn.func = "shrinkage", npn.thresh = NULL, approx = FALSE, scr = TRUE, scr.num = NULL, verbose = TRUE){	
	
	est = list()
	est$approx = approx
	
	if(is.list(L)){
		est$data = L$data
		n = nrow(L$data)
		d = ncol(L$data)
		if(!is.null(L$theta))	est$theta = L$theta
	}
	if(!is.list(L)){
		est$data = L
		n = nrow(L)
		d = ncol(L)
	}
	rm(L)
	gc(gcinfo(verbose = FALSE))	
	
	if(is.null(ind.group))	ind.group = c(1:d)
	k = length(ind.group)
		
	if(d < 3){
		cat("The fullgraph dimension < 3 and huge() will be teminated....\n")
		cat("Please refer to Pearson's product-moment correlation....\n")
		est$marker = "Terminated"
		class(est) = "huge"
		return(est)
	}
	
	if(!is.null(scr.num))
		if(scr.num == 1){
			cat("The neighborhood size < 2 and huge() will be teminated.\n")
			cat("Please refer to Pearson's product-moment correlation....\n")
			est$marker = "Terminated"
			class(est) = "huge"
			return(est)
		}
	
	
	# Nonparanormal transformation
	if(npn){
		if(is.null(npn.thresh))	npn.thresh = 1/(4*(n^0.25)*sqrt(pi*log(n)))
		est$data = huge.npn(est$data,npn.func,npn.thresh,verbose = verbose)$data
		rm(npn.thresh,npn.func)
		gc(gcinfo(verbose = FALSE))
	}
	
	if(approx){
		if(is.null(nlambda)) nlambda = 30
		if(is.null(lambda.min.ratio)) lambda.min.ratio = 0.05
		fit = huge.scr(est$data, ind.group, scr.num, approx, nlambda, lambda.min.ratio, lambda, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		rm(fit)
		gc(gcinfo(verbose = FALSE))
	}
	
	if(!approx){
		if(is.null(nlambda)) nlambda = 10
		if(is.null(lambda.min.ratio)) lambda.min.ratio = 0.1		
		if(scr){
			if(is.null(scr.num)&&(n >= d)){
				ind.mat = NULL
				est$scr = FALSE
			}
			if(is.null(scr.num)&&(n < d)){
				scr.num = min(n-1)
				ind.mat = huge.scr(est$data, ind.group = ind.group, scr.num = scr.num, verbose = verbose)$ind.mat
				est$scr = scr
			}
		}
		if(!scr){
			ind.mat = NULL
			est$scr = scr
		}
		fit = huge.subgraph(est$data, ind.group, ind.mat, alpha, lambda, nlambda, lambda.min.ratio, sym, verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		est$rss = fit$rss
		est$df = fit$df
		rm(fit)
		gc(gcinfo(verbose = FALSE))
		
		est$ind.mat = ind.mat
		est$alpha = alpha
		est$sym = sym
		rm(ind.mat,alpha,sym,approx)
		gc(gcinfo(verbose = FALSE))
	}
	est$ind.group = ind.group
	if(k<d) est$type = "Subgraph solution path"
	if(k==d) est$type = "Fullgraph solution path"
	est$npn = npn
	
	rm(ind.group,n,k,d,npn,scr,lambda,lambda.min.ratio,nlambda)
	gc(gcinfo(verbose = FALSE))	
	
	est$marker = "Successful"
	
	class(est) = "huge"
	return(est)
}

print.huge = function(x, ...){
	
	if(x$marker == "Terminated"){
		cat("huge() has been terminated\n")
		return("Please refer to the manual")
	}
	
	if(x$approx) cat("Model: Graph Approximation via Correlation Thresholding (GACT) --->",x$type,"\n")
	
	if(!x$approx){
		if(x$alpha <1) cat("Model: Meinshausen &Buhlmann Graph Estimation via Elastic Net --->",x$type,"\n")
		if(x$alpha ==1) cat("Model: Meinshausen &Buhlmann Graph Estimation via Lasso (GEL) --->",x$type,"\n")
	}
	
	if(x$npn) cat("Nonparanormal (NPN) transformation: on\n")
	
	if(!x$approx&&x$scr) cat("Graph SURE Screening (GSS): on\n")

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
	if(object$approx) cat("Model: Graph Approximation via Correlation Thresholding (GACT) --->",object$type,"\n")
	
	if(!object$approx){
		if(object$alpha <1) cat("Model: Meinshausen &Buhlmann Graph Estimation via Elastic Net --->",object$type,"\n")
		if(object$alpha ==1) cat("Model: Meinshausen &Buhlmann Graph Estimation via Lasso (GEL)--->",object$type,"\n")
	}
	
	if(object$npn) cat("Nonparanormal (NPN) transformation: on\n")
	if(!object$approx&&object$scr) cat("Graph SURE Screening (GSS): on\n")

	if(is.null(object$theta)) cat("True graph: not included\n")
	if(!is.null(object$theta)) cat("True graph: included\n")
	
	cat("Path length:",length(object$lambda),"\n")
	cat("Subgraph dimension:",length(object$ind.group),"\n")
	cat("Fullgraph dimension:",ncol(object$data),"\n")
	cat("Sparsity level:",min(object$sparsity),"----->",max(object$sparsity),"\n")
}

plot.huge = function(x, align = FALSE, ...){
	if(x$marker == "Terminated"){
		cat("huge() has been terminated\n")
		return("Please refer to the manual")
	}
	
	z.max = max(x$sparsity)
	z.min = min(x$sparsity)
	z = z.max - z.min
	z1 = which(x$sparsity>=(z.min + 0.03*z))[1]
	z2 = which(x$sparsity>=(z.min + 0.07*z))[1]
	z3 = which(x$sparsity>=(z.min + 0.15*z))[1]
	
	par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameters", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Sparsity vs. Regularization")
	
	lines(x$lambda[c(z1,z2,z3)],x$sparsity[c(z1,z2,z3)],type = "p")
	
	if(align){
		g1 = graph.adjacency(as.matrix(x$path[[z1]]), mode="undirected", diag=FALSE)
		g2 = graph.adjacency(as.matrix(x$path[[z2]]), mode="undirected", diag=FALSE)
		g3 = graph.adjacency(as.matrix(x$path[[z3]]), mode="undirected", diag=FALSE)
		layout.grid3 = layout.fruchterman.reingold(g3)
	
		plot(g1, layout=layout.grid3, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[z1],3)),sep = ""))
	
		plot(g2, layout=layout.grid3, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[z2],3)),sep = ""))
	
		plot(g3, layout=layout.grid3, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[z3],3)),sep = ""))
	}
	
	if(!align){
		g1 = graph.adjacency(as.matrix(x$path[[z1]]), mode="undirected", diag=FALSE)
		layout.grid1 = layout.fruchterman.reingold(g1)
		
		g2 = graph.adjacency(as.matrix(x$path[[z2]]), mode="undirected", diag=FALSE)
		layout.grid2 = layout.fruchterman.reingold(g2)
	
		g3 = graph.adjacency(as.matrix(x$path[[z3]]), mode="undirected", diag=FALSE)
		layout.grid3 = layout.fruchterman.reingold(g3)
	
		plot(g1, layout=layout.grid1, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[z1],3)),sep = ""))
	
		plot(g2, layout=layout.grid2, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[z2],3)),sep = ""))
	
		plot(g3, layout=layout.grid3, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[z3],3)),sep = ""))
	}
}