#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# lasso.stars(): StARS for lasso                                        #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Nov 12th, 2010                                                  #
# Verision: 0.8                                                         #
#-----------------------------------------------------------------------#
lasso.stars = function(x, y, rep.num = 20, lambda = NULL, n.lambda = 100, lambda.min = 0.001, stars.thresh = 0.1, sample.ratio = NULL, alpha = 1, verbose = TRUE){
	
	est = list()
	if(!is.matrix(x)){
		cat("Only 1 variable and lasso.stars() will be teminated....\n")
		cat("Please refer to Pearson's product-moment correlation....\n")
		est$marker = "Terminated"
		class(est) = "stars"
		return(est)
	}
	
	n = nrow(x)
	d = ncol(x)
	
	x = scale(x)
	y = y - mean(y)
	
	if(is.null(lambda)){
		lambda.max = max(abs(t(y)%*%x))/(alpha*n)
		lambda = lambda.max*exp(seq(log(1),log(lambda.min),length = n.lambda))
	}
		
	if(is.null(sample.ratio)){
		if(n>144) sample.ratio = 10*sqrt(n)/n
		if(n<=144) sample.ratio = 0.8
	}
	
	fit = glmnet(x, y, lambda = lambda)

	n.lambda = length(lambda)
	
	R.path = list()
	for(k in 1:rep.num){
    	ind.sample = sample(c(1:n), floor(n*sample.ratio), replace=FALSE)
    	out.subglm = glmnet(x[ind.sample,],y[ind.sample], lambda = fit$lambda, alpha = alpha)
    	R.path[[k]] = out.subglm$beta
    	rm(out.subglm)
		gc(gcinfo(verbose = FALSE)) 
	}
	
	if(verbose) cat('Statbility selection....')
	
	P = matrix(0,d,n.lambda)
	for(j in 1:rep.num)	P = P + abs(sign(R.path[[j]]))
	P = P/rep.num
	V = 2*apply(P*(1-P),2,sum)/d
	rm(P,rep.num)
	gc(gcinfo(verbose = FALSE))
	
	for(i in 1:n.lambda) V[i] = max(V[1:i])
	opt.index = which.min(V<=stars.thresh)
	rm(n.lambda,R.path)
	gc(gcinfo(verbose = FALSE))

	if(verbose) cat('done.\n')
	
	est$path = fit$beta
	est$lambda = lambda
	est$df = fit$df
	est$opt.index = opt.index
	est$opt.beta = fit$beta[,opt.index]
	est$opt.lambda = lambda[opt.index]
	est$Varibility = V
	est$marker = "Successful"
	class(est) = "stars"
	return(est)
}

print.stars = function(x, ...){
	if(x$marker == "Terminated"){
		cat("lasso.stars() has been terminated\n")
		return("Please refer to the manual")
	}
	cat("Path d.f.'s: from", min(x$df), "to", max(x$df),"\n")
	cat("Path length:",length(x$lambda),"\n")
}

summary.stars = function(object, ...){
	if(object$marker == "Terminated"){
		cat("lasso.stars() has been terminated\n")
		return("Please refer to the manual")
	}
	cat("Path d.f.'s: from", min(object$df), "to", max(object$df),"\n")
	cat("Path length:",length(object$lambda),"\n")
}

plot.stars = function(x, ...){
	if(x$marker == "Terminated"){
		cat("lasso.stars() has been terminated\n")
		return("Please refer to the manual")
	}
	par(mfrow = c(1,1))
	d = nrow(x$path)
	n.lambda = length(x$lambda)
	plot(x$lambda, x$path[1,], log = "x", xlab = "Regularization Parameter", ylab = "Coefficients", type = "l",xlim = rev(range(x$lambda)),ylim = c(min(x$path),max(x$path)), main = "Regularization path")
	for(j in 2:d)	lines(x$lambda, x$path[j,],xlim = rev(range(x$lambda)))
	lines(c(x$opt.lambda,x$opt.lambda),c(min(x$path),max(x$path)))
}
