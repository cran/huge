#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation (HUGE)          #
# huge(): Draw ROC Curve for a solution path                            #
#         Must have a ground truth                                      #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tourzhao@andrew.cmu.edu>; <hanliu@cs.jhu.edu>                #
# Date: Feb 28th 2010                                                   #
# Version: 1.0                                                          #
#-----------------------------------------------------------------------#

huge.roc = function(est, theta, verbose = TRUE){
	gcinfo(verbose = FALSE)
	if(class(est) == "huge"){
		G = est$path
		if(is.null(est$theta)) cat("Error, true graph is not included!\n")
		theta = est$theta
	}
	if(class(est) != "huge") G = est
	rm(est)
	gc()
	
	ROC = list()
	
	theta = as.matrix(theta)
	d = ncol(theta)
	pos.total = sum(theta!=0)
	neg.total = d*(d-1) - pos.total
	
	if(verbose) cat("Computing F1 scores, false positive rates and true positive rates....")
	ROC$tp = rep(0,length(G))
   	ROC$fp = rep(0,length(G))
   	ROC$F1 = rep(0,length(G))
   	for (r in 1:length(G)){
   		tmp = as.matrix(G[[r]]) 
   		tp.all = (theta!=0)*(tmp!=0)
   		diag(tp.all) = 0
		ROC$tp[r] <- sum(tp.all!=0)/pos.total
		fp.all = (theta==0)*(tmp!=0)
		diag(fp.all) = 0
		ROC$fp[r] <- sum(fp.all!=0)/neg.total
		
		fn = 1 - ROC$tp[r]
		precision = ROC$tp[r]/(ROC$tp[r]+ROC$fp[r])
		recall = ROC$tp[r]/(ROC$tp[r]+fn)
		ROC$F1[r] = 2*precision*recall/(precision+recall)
		if(is.na(ROC$F1[r]))	ROC$F1[r] = 0
	}
	if(verbose) cat("done.\n")
		
	rm(precision,recall,tp.all,fp.all,G,theta,fn)
   	gc()	
		
	ord.fp = order(ROC$fp)
	
	tmp1 = ROC$fp[ord.fp]
	tmp2 = ROC$tp[ord.fp]
	par(mfrow = c(1,1))	
	plot(tmp1,tmp2,type="b",main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate",ylim = c(0,1))
	ROC$AUC = sum(diff(tmp1)*(tmp2[-1]+tmp2[-length(tmp2)]))/2
	
	rm(ord.fp, tmp1, tmp2)
	gc()
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