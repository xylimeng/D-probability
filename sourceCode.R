# source code for D-probabilities 

# source function for efficient Gaussian process fitting 
# Note: 
#  - use flat prior 
#  - return type 1 and type 2 KL

library(BAS)

### Utility function for Gaussian process fitting & returning KLs for one design matrix ----
# TR calculates the tract of a squared matrix 
tr <- function(A) {sum(diag(A))}

# anisotropic: does not affect this function; keep it here for consistency of arguments between this function and its derivative 
neg.log.p.pseduo <- function(hyper, y, a, b, pre.kernel, jitter = 0, anisotropic = TRUE){
  log.r.sq = hyper[1]; log.lambda.sq = hyper[-1]; 
  n = length(y)
  lambda = exp(log.lambda.sq/2)
  r.sq = exp(log.r.sq)
  
  Q = exp(-apply(pre.kernel, c(1,2), function(k) sum(k / (lambda^2)))) * r.sq + 
    diag(1 + jitter, nrow = n, ncol = n) # jittered for stability (not quite useful if we use tau^2 = r^2 * sigma^2) 
  alpha.star = n/2 + a
  beta.star = c((y %*% solve(Q, y)))/2 + b
  # sigma.sq.hat = c((y %*% solve(Q, y)) / n) 
  log.p = -0.5 * (2 * alpha.star * log(beta.star) + determinant(Q, logarithm = TRUE)$modulus)
  return(-c(log.p))
}

# if anisotropic = TRUE: differet smoothing parameter for each covariate
# if FALSE: the same smoothing parameter 
D.neg.log.p.pseduo <- function(hyper, y, a, b, pre.kernel, jitter = 0, 
                               anisotropic = TRUE){
  log.r.sq = hyper[1]; log.lambda.sq = hyper[-1]; d = length(log.lambda.sq)
  n = length(y)
  lambda = exp(log.lambda.sq/2)
  r.sq = exp(log.r.sq)
  
  Q1 = exp(-apply(pre.kernel, c(1,2), function(k) sum(k / (lambda^2)))) * r.sq 
  Q = Q1 + diag(1 + jitter, nrow = n, ncol = n) # jittered for stability 
  
  alpha.star = n/2 + a
  beta.star = c((y %*% solve(Q, y)))/2 + b
  
  D.log.r.sq = Q1
  alpha = solve(Q, y)
  part1 = (alpha %*% t(alpha))/(beta.star/alpha.star) - solve(Q)
  
  part2 = sapply(1:d, function(k) tr(part1 %*% (Q1 * pre.kernel[,,k]))/(lambda[k]^2))
  if (!anisotropic){
    part2 = rep(sum(part2), d) 
  }
  
  D = c(tr(part1 %*% D.log.r.sq), part2)
  D = D/2 
  return(-D)
}

## GP.fit(x, y, a, b, jitter = 0, hyper.start = NULL, anisotropic = TRUE)
GP.fit <- function(x, y, a, b, jitter = 0, hyper.start = NULL, anisotropic = TRUE){
  
  x = as.matrix(x)
  
  d = ncol(x) # dimension of x exclusing the intercept 
  for (i in 1:d){
    if (var(c(x[,d])) == 0){
      x[,d] = NULL
    }
  }
  
  # pre.kernel becomes multivariate dimension 
  ## univariate
  ## pre.kernel <- outer(x, x, function(x,y) (x - y)^2/2)
  aa = matrix(0, nrow(x), nrow(x)) # template value 
  pre.kernel = vapply(1:d, function(k) outer(x[,k], x[,k], function(x, y) (x - y)^2/2), 
                      FUN.VALUE = aa)
  if (is.null(hyper.start)){
    hyper.start = c(r.sq = 0, lambda.sq = rep(0, d))
  }
  ret.pseduo <- nlminb(hyper.start, neg.log.p.pseduo, D.neg.log.p.pseduo, 
                       y = y, pre.kernel = pre.kernel, a = a, b = b, 
                       jitter = jitter, anisotropic = anisotropic, 
                       control = list(iter.max = 100, rel.tol = 8e-10, 
                                      trace = 0))
  
  hyper = ret.pseduo$par 
  log.r.sq = hyper[1]; log.lambda.sq = hyper[2:(d + 1)]; 
  n = length(y)
  lambda = exp(log.lambda.sq/2)
  r.sq = exp(log.r.sq)
  
  K = exp(-apply(pre.kernel, c(1,2), function(k) sum(k / (lambda^2)))) * r.sq
  Q = K + diag(1 + jitter, nrow = n, ncol = n) # jittered for stability 
  # sigma.sq.hat = c((y %*% solve(Q, y)))/n 
  # 'sigma.sq.hat' not useful 
  
  alpha.star = n/2 + a
  beta.star = c((y %*% solve(Q, y)))/2 + b
  sigma.sq = beta.star / (alpha.star - 1) # posterior mean 
  
  
  H = K %*% solve(Q) 
  Y.hat = H %*% y # independent of sigma.sq 
  R = sum(y * y) - sum(y * Y.hat)
  
  # for prediction 
  Q.inv.y = solve(Q, y)
  
  detail = list(H = H, Y.hat = Y.hat, alpha.star = alpha.star, 
                beta.star = beta.star, sigma.sq = sigma.sq, R = R, 
                tr.H = tr(H), Q.inv.y = Q.inv.y)
  # v = K.svd$v
  # d = K.svd$d 
  # # Y.hat = H %*% Y 
  r.sq = unname(r.sq); lambda = unname(lambda); 
  ret.par = list(tau.sq = r.sq, lambda.sq = lambda^2, detail = detail) 
  class(ret.par) <- "GP"
  return(ret.par)
} 

# x_new, x_train do NOT include intercept 
predict.GP <- function(ret_GP, x_new, x_train){
  
  ## pre.kernel <- outer(x, x, function(x,y) (x - y)^2/2)
  aa = matrix(0, nrow(x_new), nrow(x_train)) # template value 
  pre.kernel = vapply(1:ncol(x_new), 
                      function(k) outer(x_new[,k], x_train[,k], function(x, y) (x - y)^2/2), 
                      FUN.VALUE = aa)
  
  lambda.sq = ret_GP$lambda.sq
  tau.sq = ret_GP$tau.sq
  K = exp(-apply(pre.kernel, c(1,2), function(k) sum(k / lambda.sq))) * tau.sq
  
  ret = K %*% ret_GP$detail$Q.inv.y 
  
  return(ret)
}



GP.given <- function(x, y, a, b, tau.sq, lambda.sq, jitter = 0){
  
  x = as.matrix(x)
  
  d = ncol(x) # dimension of x exclusing the intercept 
  for (i in 1:d){
    if (var(c(x[,d])) == 0){
      x[,d] = NULL
    }
  }
  
  # pre.kernel becomes multivariate dimension 
  ## univariate
  ## pre.kernel <- outer(x, x, function(x,y) (x - y)^2/2)
  aa = matrix(0, nrow(x), nrow(x)) # template value 
  pre.kernel = vapply(1:d, function(k) outer(x[,k], x[,k], function(x, y) (x - y)^2/2), 
                      FUN.VALUE = aa)
  
  n = length(y)
  lambda = sqrt(lambda.sq)
  r.sq = tau.sq
  
  K = exp(-apply(pre.kernel, c(1,2), function(k) sum(k / (lambda^2)))) * r.sq
  Q = K + diag(1 + jitter, nrow = n, ncol = n) # jittered for stability 
  # sigma.sq.hat = c((y %*% solve(Q, y)))/n 
  # 'sigma.sq.hat' not useful 
  
  alpha.star = n/2 + a
  beta.star = c((y %*% solve(Q, y)))/2 + b
  sigma.sq = beta.star / (alpha.star - 1) # posterior mean 
  
  
  H = K %*% solve(Q) 
  Y.hat = H %*% y # independent of sigma.sq 
  R = sum(y * y) - sum(y * Y.hat)
  
  detail = list(H = H, Y.hat = Y.hat, alpha.star = alpha.star, 
                beta.star = beta.star, sigma.sq = sigma.sq, R = R, 
                tr.H = tr(H))
  # v = K.svd$v
  # d = K.svd$d 
  # # Y.hat = H %*% Y 
  r.sq = unname(r.sq); lambda = unname(lambda); 
  ret.par = list(tau.sq = r.sq, lambda.sq = lambda^2, detail = detail) 
  class(ret.par) <- "GP"
  return(ret.par)
} 
print.GP <- function(ret.par){
  summary = c(tau.sq = ret.par$tau.sq, 
              lambda.sq = ret.par$lambda.sq, 
              sigma.sq = ret.par$detail$sigma.sq)
  
  print(summary)
}

data2KL <- function(X.j, y, GP, g, type = 1, a = 0, b = 0){
  n = length(y); Y = y
  # H.j = X %*% solve(A) %*% t(X); used to be 'part'
  H.j = X.j %*% solve(t(X.j) %*% X.j, t(X.j))
  
  if (is.infinite(g)){
    H.j = H.j
  } else {
    H.j = g/(g + 1) * H.j 
  }
  
  Y.hat.j = H.j %*% Y
  R.j = sum(Y * Y) - sum(Y * Y.hat.j) # Y^T (I - H.j) Y 
  
  Y.hat = GP$detail$Y.hat 
  H = GP$detail$H
  tr.H = tr(GP$detail$H)
  
  diff.Y.hat = Y.hat.j - Y.hat
  alpha.star.j = n/2 + a
  beta.star.j = R.j/2 + b
  alpha.star = GP$detail$alpha.star
  beta.star = GP$detail$beta.star 
  
  # three constants: (1/sigma.j^2, sigma^2/sigma.j^2, log(sigma.j^2/sigma^2))) 
  # either use point estimates or use integral  
  const = rep(NA, 3)
  
    const[1] = alpha.star.j / beta.star.j 
    const[2] = const[1] * beta.star / (alpha.star - 1)
    const[3] = log(beta.star.j) - digamma(alpha.star.j) - 
      (log(beta.star) - digamma(alpha.star))
  
  if (type == 1){ # use posterior mean's 
    double.KL.j = mean(diff.Y.hat ^2) * const[1] + 
      (1 + tr.H / n) * const[2] + const[3] + tr(H.j)/n - 1
  }
  
  if (type == 2){ # use predictive densities conditional on (sigma_j, sigma)
    double.KL.j = mean(diff.Y.hat * solve(H.j + diag(n), diff.Y.hat)) * const[1] + 
      tr(solve(H.j + diag(n), H + diag(n)))/n * const[2] + const[3] - 1 + (determinant(H.j + diag(n))$modulus - determinant(H + diag(n))$modulus)/n
  }
  
  KL.j = double.KL.j / 2
  # calculate n * KL 
  return(KL.j)
}

all.in.one <- function(x, y, output = 1){
  GP <- GP.fit(x, y, a = 0, b = 0)
  n = length(y)
  X = cbind(1, x)
  ret = matrix(NA, 2, 2); col.name = rep(NA, 2)
  g.g = Inf # use flat priors 
  
  for (g.type in 1:2){
      idx = sprintf('KL%d', g.type)
      col.name[g.type] = idx
      ret[1, g.type] = data2KL(X, y, GP, g = g.g, type = g.type)
      ret[2, g.type] = data2KL(X[,1], y, GP,  g = g.g, type = g.type)
    
  }
  
  colnames(ret) = col.name
  
  if (output == 1){
    return_value = list(GP = GP, KL = ret, n = n)
  }
  
  if (output == 0){
    return_value = list(KL = ret, n = n)
  }
  
  return(return_value)
}


### Utility functions for many p (covariates) ------
# including both BAS and D-Bayes 

# Utility functions 
# rank starts from 1; intercept is always included 
binary2rank <- function(x){
  # x: (intercept, x1, x2, ..., xp)
  p = length(x) - 1
  sum(x[-1] * 2^(0:(p - 1))) + 1 
}

p2model.id <- function(p){
  m <- t(sapply(0:(2^p - 1),function(x){ as.integer(intToBits(x))[1:p]}))
  idx = lapply(1:2^p, function(k) which(m[k, ] == 1))
  rank = sapply(1:2^p, function(k) binary2rank(c(1, m[k, ]))) 
  
  return(list(m = m, idx = idx, rank = rank))
}

# return (-n * kl): DO NOT take the exp at this moment 
conv.summary.KL <- function(input, model.id, n, num.model){
  # num.model = min(length(input), num.model)
  idx = order(input, decreasing = FALSE)[1:num.model]
  ret = c(model.id$rank[idx], -n * input[idx])
  return(ret)
}

# input: ret = BAS.lm from method "hyper.g" 
# output: BAS from other method 
conv.update <- function(method, ret){
  if ((method == "g-prior") | (method == "hyper-g-laplace")){
    ret = update(ret, newprior=method, alpha = 3)
  } else {
    ret = update(ret, newprior=method)
  }
  ret 
}

# summarize each method 
conv.summary <- function(ret, num.model){
  p = ncol(ret$X) - 1
  rank = apply(summary(ret, n.models = num.model)[, 1:(p + 1)], 1, binary2rank)
  posterior.prob <- summary(ret, n.models = num.model)[, p + 3]
  return(c(rank, posterior.prob))
}


# add num.model 
all.in.one_D_Bayes <- function(x, y, num.model, g.a, g.b, g.g, anisotropic = TRUE){
  
  p = ifelse(is.null(ncol(x)), 1, ncol(x))
  y = as.vector(y)
  n = length(y)
  model.id = p2model.id(p)
  
  ## D-Bayes
  # for (g.ab in c(0, 0.1, 1, 10, 100)){
  ret = GP.fit(x, y, a = g.a, b = g.b, jitter = 0, anisotropic = anisotropic)
  #  mean((y - ret$detail$Y.hat)^2)
  
  KL.all = matrix(NA, nrow = 2, ncol = 2^p) # all 2^p models 
  for (g.type in 1:2){
      KL.all[g.type, ] = unlist(lapply(model.id$idx, function(k)
        data2KL(X.j = cbind(1, x[, k]), y, g = g.g , GP = ret, type = g.type,
                sigma.method = g.sigma.method, a = 0, b = 0)))
  }
  
  summary.D_Bayes = sapply(1:2, function(k) 
    conv.summary.KL(KL.all[k, ], model.id, n, num.model))
  colnames(summary.D_Bayes) = c("KL1", "KL2")
  # }
  
  ## BAS 
  all.method = c("hyper-g", "AIC", "BIC", "g-prior", "ZS-null", "ZS-full", "hyper-g-laplace", "hyper-g-n", "EB-local", "EB-global")
  
  bas_all = as.list(1:length(all.method))
  names(bas_all) <- all.method
  
  ret = bas.lm(y ~ x, method = "BAS", 
               prior = "hyper-g", alpha = 3, modelprior = uniform())
  bas_all[["hyper-g"]] <- ret
  for (method in all.method[2:length(all.method)]){
    bas_all[[method]] <- tryCatch(conv.update(method, ret), error = function(e) NA)
  }
  
  bas_summary = matrix(NA, num.model * 2, ncol = length(all.method))
  colnames(bas_summary) = all.method
  
  for (method in all.method){
    bas_summary[ ,method] = tryCatch(
      conv.summary(bas_all[[method]], num.model), error = function(e) rep(NA, 10))
  }
  
  all_summary = cbind(summary.D_Bayes, bas_summary) 
  return(all_summary)
  
}

