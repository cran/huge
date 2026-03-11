# Internal helper functions (not exported)

.huge_preprocess = function(x, verbose = TRUE) {
  n = nrow(x)
  d = ncol(x)
  cov.input = isSymmetric(x)
  if(cov.input) {
    if(verbose) cat("The input is identified as the covariance matrix.\n")
    S = cov2cor(x)
  } else {
    x = scale(x)
    S = cor(x)
  }
  list(x = x, S = S, n = n, d = d, cov.input = cov.input)
}

.huge_default_lambda = function(S, d, nlambda = NULL, lambda.min.ratio = NULL, lambda = NULL) {
  if(!is.null(lambda)) return(list(lambda = lambda, nlambda = length(lambda)))
  if(is.null(nlambda)) nlambda = 10
  if(is.null(lambda.min.ratio)) lambda.min.ratio = 0.1
  lambda.max = max(max(S - diag(d)), -min(S - diag(d)))
  lambda.min = lambda.min.ratio * lambda.max
  lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  list(lambda = lambda, nlambda = nlambda)
}
