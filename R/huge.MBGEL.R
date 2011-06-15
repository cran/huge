#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge.mbgel(): Meinshausen & Buhlmann Graph                            #
#				Estimation via Lasso (MBGEL)                            #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@gmail.com>; <hanliu@cs.jhu.edu>                     #
# Date: Jun 8th 2011                                                    #
# Version: 1.0.2                                                        #
#-----------------------------------------------------------------------#

## Main function
huge.mbgel = function(x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, scr = NULL, scr.num = NULL, idx.mat = NULL, sym = "or", verbose = TRUE)
{
	gcinfo(FALSE)
	n = nrow(x);
	d = ncol(x);
	fit = list()
	fit$cov.input = isSymmetric(x);
	if(fit$cov.input)
	{
		if(verbose) cat("The input is identified as the covriance matrix.\n")
		S = cov2cor(x);
	}
	if(!fit$cov.input)
	{
		x = scale(x)
		S = t(x)%*%x/n
	}
	
	rm(x)
	gc()
	
	if(is.null(idx.mat))
	{
		if(is.null(scr)&&(n>=d))
			scr = FALSE
		if(is.null(scr)&&(n<d))
			scr = TRUE
		
		if(scr)
		{
			if(is.null(scr.num))
			{
				if(n<d)
					scr.num = n-1
				if(n>=d)
				{
					if(verbose) cat("Graph SURE Screening (GSS) is skipped without specifying scr.num.\n")
					scr = FALSE
				}	
			}
		}
		fit$scr = scr
	}
	
	if(!is.null(idx.mat))
	{
		scr = TRUE
		fit$scr = scr
		scr.num = nrow(idx.mat)
	}	

	if(!is.null(lambda)) nlambda = length(lambda)
	if(is.null(lambda))
	{
		if(is.null(nlambda))
			nlambda = 10
		if(is.null(lambda.min.ratio))
			lambda.min.ratio = 0.1
		lambda.max = max(max(S-diag(d)),-min(S-diag(d)))
		lambda.min = lambda.min.ratio*lambda.max
		lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
		rm(lambda.max,lambda.min,lambda.min.ratio)
		gc()
	}
	maxdf = min(d,n);
   	if(scr)
   	{
   		
   		if(verbose)
   		{
   			cat("Conducting Meinshausen & Buhlmann Graph Estimation via Lasso (MBGEL) with Graph Sure Screening....")
        	flush.console()
   		}

		if(is.null(idx.mat))
			idx.mat = apply(-abs(S),2,order)[2:(scr.num+1),] - 1
		
		fit$idx.mat = idx.mat
		out=.C("SPMBscr", S=as.double(S), idx_scr=as.integer(idx.mat), lambda=as.double(lambda), nnlambda = as.integer(nlambda), dd = as.integer(d), nnscr = as.integer(scr.num), x = as.double(rep(0,d*maxdf*nlambda)),col_cnz = as.integer(rep(0,d+1)), row_idx = as.integer(rep(0,d*maxdf*nlambda)),PACKAGE="huge")
		for(i in 1:d)
		{
			if(out$col_cnz[i+1]>out$col_cnz[i])
			{
				idx.tmp = (out$col_cnz[i]+1):out$col_cnz[i+1]
				ord = order(out$row_idx[idx.tmp])
				out$row_idx[idx.tmp] = out$row_idx[ord + out$col_cnz[i]]
				out$x[idx.tmp] = out$x[ord + out$col_cnz[i]]
			}
		}
	}
   	   	
   	if(!scr)
   	{
   		if(verbose)
   		{
   			cat("Conducting Meinshausen & Buhlmann Graph Estimation via Lasso (MBGEL)....")
        	flush.console()
   		}
		fit$idx_mat = NULL
		out=.C("SPMBgraph", S=as.double(S), lambda=as.double(lambda), nnlambda = as.integer(nlambda), dd = as.integer(d), x = as.double(rep(0,d*maxdf*nlambda)),col_cnz = as.integer(rep(0,d+1)), row_idx = as.integer(rep(0,d*maxdf*nlambda)),PACKAGE="huge")
		for(i in 1:d)
		{
			if(out$col_cnz[i+1]>out$col_cnz[i])
			{
				idx.tmp = (out$col_cnz[i]+1):out$col_cnz[i+1]
				ord = order(out$row_idx[idx.tmp])
				out$row_idx[idx.tmp] = out$row_idx[ord + out$col_cnz[i]]
				out$x[idx.tmp] = out$x[ord + out$col_cnz[i]]
			}
		}
	}
	
	   	
   	G = new("dgCMatrix", Dim = as.integer(c(d*nlambda,d)), x = as.vector(out$x[1:out$col_cnz[d+1]]),p = as.integer(out$col_cnz), i = as.integer(out$row_idx[1:out$col_cnz[d+1]]))

	fit$beta = list()
	fit$path = list()
	fit$df = matrix(0,d,nlambda)
	fit$rss = matrix(0,d,nlambda)	
	fit$sparsity = rep(0,nlambda)	
	for(i in 1:nlambda)
	{
		fit$beta[[i]] = G[((i-1)*d+1):(i*d),]
		fit$path[[i]] = abs(fit$beta[[i]])
		fit$df[,i] = apply(sign(fit$path[[i]]),2,sum)
		
		if(sym == "or")
			fit$path[[i]] = sign(fit$path[[i]] + t(fit$path[[i]]))
		if(sym == "and")
			fit$path[[i]] = sign(fit$path[[i]] * t(fit$path[[i]]))
		fit$sparsity[i] = sum(fit$path[[i]])/d/(d-1)
	}	
   	rm(G,idx.tmp,ord,out)
   	
   	fit$lambda = lambda

	if(verbose)
   	{
   		cat("done\n")
        flush.console()
    }
   		
   	rm(verbose,nlambda)
   	gc()
   	return(fit)
}