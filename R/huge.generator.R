#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.generator(): The data generating function                        #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@gmail.com>; <hanliu@cs.jhu.edu>                     #
# Date: Nov 9th, 2010                                                   #
# Version: 0.7
#-----------------------------------------------------------------------#
huge.generator = function(n = 200, d = 50, graph = "random", v = NULL, u = NULL, g = NULL, prob = 0.03, vis = FALSE, verbose = TRUE){	
	if(verbose) cat("Generating data from the multivariate normal distribution with the", graph,"graph structure....")
	if(is.null(g)){
		g = 1
		if(graph == "hub"){
			if(d >= 16)	g = floor(d/8)
			if(d < 16) g = 2
		}
		if(graph == "clique"){
			if(d >= 10)	g = floor(d/5)
			if(d < 10) g = 2
		}
	}	
	# parition variables into groups
	g.large = d%%g
	g.small = g - g.large
	n.small = floor(d/g)
	n.large = n.small+1
	g.list = c(rep(n.small,g.small),rep(n.large,g.large))
	g.ind = rep(c(1:g),g.list)
	
	rm(g.large,g.small,n.small,n.large,g.list)
	gc(gcinfo(verbose = FALSE))
	
	# build the graph structure
	theta = matrix(0,d,d);
	if(graph == "band"){
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.1
		for(i in 1:g){
			diag(theta[1:(d-i),(1+i):d]) = 1
			diag(theta[(1+i):d,1:(d-1)]) = 1
		}	
	}
	if(graph == "clique"){
		if(is.null(u)) u = 0
		if(is.null(v)) v = 0.3
		for(i in 1:g){
		 	tmp = which(g.ind==i)
		 	theta[tmp,tmp] = 1
		 	rm(tmp)
		 	gc(gcinfo(verbose = FALSE))
		}
	}
	if(graph == "hub"){
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.3
		for(i in 1:g){
		 	tmp = which(g.ind==i)
		 	theta[tmp[1],tmp] = 1
		 	theta[tmp,tmp[1]] = 1
		 	rm(tmp)
		 	gc(gcinfo(verbose = FALSE))
		}
	}
	if(graph == "random"){
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.3
		tprob = sqrt(prob/2)*(prob<0.5) + (1-sqrt(0.5-0.5*prob))*(prob>=0.5) 
		tmp = matrix(runif(d^2,0,0.5),d,d)
		tmp = tmp + t(tmp)
		theta[tmp < tprob] = 1
		theta[tmp >= tprob] = 0
		rm(tmp)
		gc(gcinfo(verbose = FALSE))
	}
	diag(theta) = 0
	omega = theta*v
	
	# make omega positive definite and standardized
	diag(omega) = abs(min(eigen(omega)$values)) + 0.01 + u
	sigma = cov2cor(solve(omega))
	omega = solve(sigma)
	
	# generate multivariate normal data
	x = mvrnorm(n,rep(0,d),sigma)
	#x = scale(x)
	
	# graph and covariance visulization
	if(vis == TRUE){
		fullfig = par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
		fullfig[1] = image(theta, col = gray.colors(256),  main = "Adjacency Matrix")

		fullfig[2] = image(sigma, col = gray.colors(256), main = "Covariance Matrix")
		g = graph.adjacency(theta, mode="undirected", diag=FALSE)
		layout.grid = layout.fruchterman.reingold(g)

		fullfig[3] = plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA,main = "Graph Pattern")

		fullfig[4] = image(omega, col = gray.colors(256), main = "Precision Matrix")
		rm(fullfig,g,layout.grid)
		gc(gcinfo(verbose = FALSE))
	}
	if(verbose) cat("done.\n")
	rm(vis,verbose)
	gc(gcinfo(verbose = FALSE))
	
	sim = list(data = x, sigma = sigma, omega = omega, theta = Matrix(theta,sparse = TRUE), n=n, d=d, sparsity= sum(theta)/(d*(d-1)), graph.type=graph)
	class(sim) = "sim" 
	return(sim)
}


print.sim = function(x, ...){
	cat("The dimension of the generated data matrix: n =", x$n, "d=", x$d, "\n")
   cat("Graph type = ", x$graph.type, "\n")
}

summary.sim = function(object, ...){
	cat("The dimension of the generated data matrix: n =", object$n, "d=", object$d, "\n")
	cat("Graph type = ", object$graph.type, "\n")
}

plot.sim = function(x, ...){	
   oldpar = par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
   image(as.matrix(x$theta), col = gray.colors(256),  main = "Adjacency Matrix")
	image(x$sigma, col = gray.colors(256), main = "Covariance Matrix")
	g = graph.adjacency(x$theta, mode="undirected", diag=FALSE)
	layout.grid = layout.fruchterman.reingold(g)
	
	plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=5, vertex.label=NA,main = "Graph Pattern")
	image(x$omega, col = gray.colors(256), main = "Precision Matrix")
	rm(g, layout.grid)
	gc(gcinfo(verbose = FALSE))
	par(oldpar)
}

