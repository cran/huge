#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge(): Draw ROC Curve for a solution path                            #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@gmail.com>; <hanliu@cs.jhu.edu>                     #
# Date: Nov 9th 2010                                                    #
# Version: 0.7                                                          #
#-----------------------------------------------------------------------#

huge.roc = function(est, theta, ind.group, verbose = TRUE){
	
	if(class(est) == "huge"){
		G = est$path
		ind.group = est$ind.group
		if(is.null(est$theta)) cat("Error, true graph is not included!\n")
		theta = est$theta
	}
	if(class(est) != "huge") G = est
	rm(est)
	gc(gcinfo(verbose = FALSE))
	
	ROC = list()
	
	if(missing(ind.group)) ind.group = c(1:ncol(theta))
	k = length(ind.group)
	
	sub.theta = as.matrix(theta[ind.group,ind.group])
	rm(theta)
   	gc(gcinfo(verbose = FALSE))	
	
	sub.G = list()
	for(i in 1:length(G)) sub.G[[i]] = G[[i]][ind.group,ind.group]
	rm(G,ind.group)
   	gc(gcinfo(verbose = FALSE))
	
	pos.total = sum(sub.theta!=0)
	neg.total = k*(k-1) - pos.total
	
	if(verbose) cat("Computing F1 scores, false positive rates and true positive rates....")
	ROC$tp = rep(0,length(sub.G))
   	ROC$fp = rep(0,length(sub.G))
   	ROC$F1 = rep(0,length(sub.G))
   	for (r in 1:length(sub.G)){
   		tmp = as.matrix(sub.G[[r]]) 
   		tp.all = (sub.theta!=0)*(tmp!=0)
   		diag(tp.all) = 0
		ROC$tp[r] <- sum(tp.all!=0)/pos.total
		fp.all = (sub.theta==0)*(tmp!=0)
		diag(fp.all) = 0
		ROC$fp[r] <- sum(fp.all!=0)/neg.total
		
		fn = 1 - ROC$tp[r]
		precision = ROC$tp[r]/(ROC$tp[r]+ROC$fp[r])
		recall = ROC$tp[r]/(ROC$tp[r]+fn)
		ROC$F1[r] = 2*precision*recall/(precision+recall)
		if(is.na(ROC$F1[r]))	ROC$F1[r] = 0
	}
	if(verbose) cat("done.\n")
		
	rm(precision,recall,tp.all,fp.all,sub.G,sub.theta,fn)
   	gc(gcinfo(verbose = FALSE))	
		
	ord.fp = order(ROC$fp)
	
	tmp1 = ROC$fp[ord.fp]
	tmp2 = ROC$tp[ord.fp]
	par(mfrow = c(1,1))	
	plot(tmp1,tmp2,type="b",main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate",ylim = c(0,1))
	ROC$AUC = sum(diff(tmp1)*(tmp2[-1]+tmp2[-length(tmp2)]))/2
	
	rm(ord.fp, tmp1, tmp2)
	gc(gcinfo(verbose = FALSE))
	class(ROC) = "roc"
	return(ROC)
}

print.roc = function(x, ...){
	cat("True Postive Rate: from",min(x$tp),"to",max(x$tp),"\n")
	cat("False Positive Rate: from",min(x$fp),"to",max(x$fp),"\n")
	cat("Area under Curve:",x$AUC,"\n")
	cat("Maximum F1 Score:",max(x$F1),"\n")
}

summary.roc = function(object, ...){
	cat("True Postive Rate: from",min(object$tp),"to",max(object$tp),"\n")
	cat("False Positive Rate: from",min(object$fp),"to",max(object$fp),"\n")
	cat("Area under Curve:",object$AUC,"\n")
	cat("Maximum F1 Score:",max(object$F1),"\n")
}

plot.roc = function(x, ...){	
	ord.fp = order(x$fp)
	par(mfrow = c(1,1))
	plot(x$fp[ord.fp],x$tp[ord.fp],type="b",main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate",ylim = c(0,1))
}