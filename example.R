# Example 

rm(list = ls())
source("sourceCode.R")

### UScrime data 
library(MASS)
data("UScrime")

y = log(UScrime[, 16])
x = as.matrix(UScrime[, 1:15])
x[, -2] = log(x[, -2])

# rescale x: mean 0, sd 1
rescale_x = x
for (i in 1:ncol(x)){
    rescale_x[,i] = (x[, i] - min(x[,i])) / (max(x[,i]) - min(x[,i]))
}
x = rescale_x

# fit the reference model 
g.ab = 0 
ret = GP.fit(x, y, a = g.ab, b = g.ab, jitter = 0)
# mean((y - ret$detail$Y.hat)^2) # residual sum of squares 

n = length(y)
p = ncol(x)
m <- t(sapply(1:(2^p - 1),function(x){ as.integer(intToBits(x))[1:p]}))
idx = lapply(1:(2^p - 1), function(k) which(m[k, ] == 1))

KL.all = matrix(NA, nrow = 4, ncol = 2^p) # null model is the last column
count = 0
for (g.type in 1:2){
  KL.all[g.type, 1:(2^p - 1)] = unlist(lapply(idx, function(k)
    data2KL(X.j = cbind(1, x[, k]), y, g = Inf, GP = ret, type = g.type, a = g.ab, b = g.ab)))
}

# add the null model 



### UScrime data ------
repar <- function(prior.mean, prior.sd){
  prior.var = prior.sd^2
  g.a = prior.mean^2 / prior.var + 2
  g.b = prior.mean * (g.a - 1)
  
  # sprintf('prior mean = (%.4f); 
  #         prior standard deviation = (%.4f)', 
  #         g.b / (g.a - 1), sqrt(g.b^2 / (g.a - 1)^2 / (g.a - 2)))
  
  c(shape = g.a, rate = g.b)
}

y = dat$uscrime$y; x = dat$uscrime$x.rescale
lm_fit = lm(y ~ x)
prior.var = summary(lm_fit)$sigma^2

g.ab = repar(prior.var, 0.02)

p = ncol(x)
ret_uscrime_1 = all.in.one_D_Bayes(x = x, y = y, num.model = 2^p, g.a = g.ab[1], 
                                   g.b = g.ab[2], g.g = Inf, anisotropic = FALSE)

ret_uscrime_2 = all.in.one_D_Bayes(x = x, y = y, num.model = 2^p, g.a = g.ab[1], 
                                   g.b = g.ab[2], g.g = length(y), anisotropic = FALSE)

y = dat$uscrime$y; x = dat$uscrime$x
lm_fit = lm(y ~ x)
prior.var = summary(lm_fit)$sigma^2

g.ab = repar(prior.var, 0.02)

p = ncol(x)
ret_uscrime_3 = all.in.one_D_Bayes(x = x, y = y, num.model = 2^p, g.a = g.ab[1], 
                                   g.b = g.ab[2], g.g = Inf, anisotropic = FALSE)

ret_uscrime_4 = all.in.one_D_Bayes(x = x, y = y, num.model = 2^p, g.a = g.ab[1], 
                                   g.b = g.ab[2], g.g = length(y), anisotropic = FALSE)




#################################################
##############Summary: ozone #####################
#################################################
# g-prior: g = n 
# hyper-g-prior: Liang et al. 
rm(list = ls())
load("ret.RData")

my.plot <- function(x, y, ...){
  upper = max(max(x), max(y))
  plot(x, y, xlim = c(0, upper), ylim = c(0, upper), ...)
  abline(c(0, 1), lwd = 2)
}

## Effect of prior specification
## g = n gives ridicularly small D-prob. for small sample size 
## because everthing is in the exponential 
## For large sample size, g doesn't make any difference 
## # weird problem: Kl1s, KL2s, KL3s - all 0 - solved by rescaling 

x = dat$ozone$x.rescale
y = dat$ozone$y

p = ncol(x)
num.model = 2^p
method.name <- c("KL1", "KL1s", "KL2", "KL2s", "KL3", "KL3s", "g-prior", "hyper-g")

Z_all = list(1:2) # two priors where g = Inf vs g = n

ret = ret_ozone_1; 
head(ret[(num.model + 1):(2 * num.model), ])
Z = matrix(NA, nrow = num.model, ncol = length(method.name))
colnames(Z) = method.name

for (i in method.name){
  model.id = as.integer(ret[1:num.model, i])
  model.prob = ret[(num.model + 1):(2 * num.model), i]
  Z[model.id, i] <- model.prob
}
Z[, 1:6] = exp(Z[, 1:6])
Z_all[[1]] = Z

ret = ret_ozone_2; 
head(ret[num.model:(2 * num.model), ])
Z = matrix(NA, nrow = num.model, ncol = length(method.name))
colnames(Z) = method.name

for (i in method.name){
  model.id = as.integer(ret[1:num.model, i])
  model.prob = ret[(num.model + 1):(2 * num.model), i]
  Z[model.id, i] <- model.prob
}
Z[, 1:6] = exp(Z[, 1:6])
Z_all[[2]] = Z


my.plot(Z_all[[1]][, "KL2s"], Z_all[[2]][, "KL2s"])


# ozone data - KL1 vs. (KL2, g-prior, hyper-g prior, Kl1 with g = n)
# all use conditional 
# 
conditional_Z = t(t(Z_all[[1]][, 1:6]) / apply(Z_all[[1]][, 1:6], 2, sum))

pdf("KL1.vs.others.pdf")
my.plot <- function(x, y, ...){
  upper = max(max(x), max(y))
  plot(x, y, xlim = c(0, upper), ylim = c(0, upper), ...)
  abline(c(0, 1), lwd = 2)
}
# plot 1: g-prior; 2: hyper-g prior; 3: KL2; 4: KL1 with g = n 
my.plot(conditional_Z[, "KL1"], Z_all[[1]][, "g-prior"], pch = 16, 
        cex.axis = 2.5, lwd = 3, xlab = "Conditional D-bayes", ylab = "g prior")
my.plot(conditional_Z[-29, "KL1"], Z_all[[1]][-29, "hyper-g"], pch = 16, 
        cex.axis = 2.5, lwd = 3, xlab = "Conditional D-bayes", ylab = "g prior")
my.plot(conditional_Z[, "KL1"], conditional_Z[, "KL2"], pch = 16, 
        cex.axis = 2, lwd = 3, xlab = "Type 1", ylab = "Type 2")
my.plot(Z_all[[1]][, "KL1"] / sum(Z_all[[1]][, "KL1"]), Z_all[[2]][, "KL1"] / sum(Z_all[[2]][, "KL1"]), cex.axis = 2, lwd = 3, pch = 16, 
        xlab = "g = infinity", ylab = "g = n")
dev.off()

pdf("KL1.vs.others.xyswitch.pdf") # switch x and y 
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 2, 0, 0)) 
my.plot <- function(x, y, xlab, ylab, ...){
  upper = max(max(x), max(y))
  plot(y, x, xlab = ylab, ylab = xlab, xlim = c(0, upper), ylim = c(0, upper), 
       pch = 19, mgp=c(6.5,1.5,0), ...)
  abline(c(0, 1), lwd = 3)
}
# par(mar=c(5,4,4,4)+0.1)
this.cex = 2.5; 
# plot 1: g-prior; 2: hyper-g prior; 3: KL2; 4: KL1 with g = n 
my.plot(conditional_Z[, "KL1"], Z_all[[1]][, "g-prior"], 
        cex.axis = this.cex, lwd = 3, xlab = "Conditional D-bayes", ylab = "g prior", las = 1)
my.plot(conditional_Z[, "KL1"], Z_all[[1]][, "hyper-g"], 
        cex.axis = this.cex, lwd = 3, xlab = "Conditional D-bayes", ylab = "g prior", las = 1)
my.plot(conditional_Z[, "KL1"], conditional_Z[, "KL2"], 
        cex.axis = this.cex, lwd = 3, xlab = "Type 1", ylab = "Type 2", las = 1)
my.plot(Z_all[[1]][, "KL1"] / sum(Z_all[[1]][, "KL1"]), Z_all[[2]][, "KL1"] / sum(Z_all[[2]][, "KL1"]), cex.axis = this.cex, lwd = 3, 
        xlab = "g = infinity", ylab = "g = n", las = 1)
rm(my.plot)
dev.off()

## which model is selected? 



my.plot(Z_all[[1]][, "KL1"], Z_all[[1]][, "g-prior"], pch = 16, 
        cex.axis = 2.5, lwd = 3, xlab = "D-bayes", ylab = "g prior")




my.plot(conditional_Z[, "KL1"], conditional_Z[, "KL1s"], pch = 16, 
        cex.axis = 2.5, lwd = 3, xlab = "Type 1", ylab = "Type 1 - simplified")


# ozone data - effect of priors 
par(mfrow = c(1, 2))
my.plot(Z_all[[1]][, "g-prior"], Z_all[[2]][, "hyper-g"], pch = 16, 
        cex.axis = 2.5, lwd = 3, xlab = "g-prior", ylab = "hyper-g prior")

my.plot(Z_all[[1]][, "KL1"], Z_all[[2]][, "KL1"], cex.axis = 2.5, lwd = 3, pch = 16, 
        xlab = "g = infinity", ylab = "g = n")

# paper plot replicates 
par(mfrow = c(1, 1))
pdf('prior.sensitivity.g.pdf')
my.plot(Z_all[[1]][, "g-prior"], Z_all[[2]][, "hyper-g"], pch = 16, 
        cex.axis = 2, lwd = 3, xlab = "g-prior", ylab = "hyper-g prior")
dev.off()

pdf('prior.sensitivity.abs.D.bayes.pdf')
# abs D-bayes 
my.plot(Z_all[[1]][, "KL1"], Z_all[[2]][, "KL1"], cex.axis = 2, lwd = 3, pch = 16,
        xlab = "g = infinity", ylab = "g = n")
dev.off()

model.id = p2model.id(p)
model.id$idx[[which.max(Z_all[[1]][, "KL1"])]]
model.id$idx[[which.max(Z_all[[1]][, "KL2"])]]
model.id$idx[[which.max(Z_all[[1]][, "g-prior"])]]
model.id$idx[[which.max(Z_all[[1]][, "hyper-g"])]]

# combine four posterior probabilities 
post.prob_all = cbind(conditional_Z[, "KL1"], conditional_Z[, "KL2"], 
                      Z_all[[1]][, "g-prior"], Z_all[[1]][, "hyper-g"])
# table: 5 by 3
# each row is for a method (four methods + reference mode)
# columns: selected model; model probability; RMSE (highest model probility)

table.str <- c("$\\pi_{j, 1 \\mid \\MM}$", "$\\pi_{j, 2 \\mid \\MM}$", "Zellner-Siow $g$ prior", "hyper-$g$ prior")
for (ith.method in 1:4){
  ith.str = table.str[ith.method]
  post.prob = post.prob_all[, ith.method]
  flag_highest = model.id$idx[[which.max(post.prob)]]
  
  output = sprintf('%s & %s & %.2f \n', 
                   ith.str, 
                   paste(colnames(x)[flag_highest], collapse = ','), 
                   post.prob[which.max(post.prob)]
  )
  table.str[ith.method] = output
}

cat(table.str)


# flag_median = which(apply(model.id$m * post.prob, 2, sum) > 0.5)
flag_highest = model.id$idx[[which.max(post.prob)]]
# flag_median
flag_highest

### compare MSE & RISK  ----
# # three methods: bas.lm with hyper-g, bas.lm with g-rpiro, and GP 
# # five terminal methods depending on "MPM" and "HPM" 
# 
# t1 = proc.time()
# N.rep = 100
# 
# MSE_all = matrix(NA, nrow = N.rep, ncol = 7)
# RISK_all = matrix(NA, nrow = N.rep, ncol = 7)
# ret_all = as.list(1:N.rep)
# data_all = as.list(1:N.rep)
# 
# ozone = data.frame(x, y)
# summarize <- function(object, ozone, x_train_id){
# 
#   if (class(object) == "bas"){
#     MSE = rep(NA, 3)
#     RISK = rep(NA, 3)
#     yhat1 = predict(object, newdata = ozone[-x_train_id, ], estimator = "MPM")$Ybma
#     yhat2 = predict(object, newdata = ozone[-x_train_id, ], estimator = "HPM")$Ybma
#     yhat3 = predict(object, newdata = ozone[-x_train_id, ], estimator = "BMA")$Ybma
#     MSE[1] = mean((ozone[-x_train_id, 'y'] - yhat1)^2)
#     MSE[2] = mean((ozone[-x_train_id, 'y'] - yhat2)^2)
#     MSE[3] = mean((ozone[-x_train_id, 'y'] - yhat3)^2)
# 
#     yhat1 = predict(object, estimator = "MPM")$fit
#     yhat2 = predict(object, estimator = "HPM")$fit
#     yhat3 = predict(object, estimator = "BMA")$fit
#     RISK[1] = mean((ozone[x_train_id, 'y'] - yhat1)^2)
#     RISK[2] = mean((ozone[x_train_id, 'y'] - yhat2)^2)
#     RISK[3] = mean((ozone[x_train_id, 'y'] - yhat3)^2)
#   }
# 
#   if (class(object) == "GP"){
#     yhat_GP = predict(object, x_train = x[x_train_id, 1:8],
#                       x_new = x[-x_train_id, 1:8])
#     MSE = mean((ozone[-x_train_id, 'y'] - yhat_GP)^2)
# 
#     yfit_GP = object$detail$Y.hat
#     RISK = mean((ozone[x_train_id, 'y'] - yfit_GP)^2)
#   }
# 
#   return(list(MSE = MSE, RISK = RISK))
# }
# 
# set.seed(20161)
# for (ith.rep in 1:N.rep)
# {
#   n = length(y)
#   x_train_id = sample(1:n, round(0.5 * n))
# 
#   ret1 = bas.lm(y ~ vh + wind + humidity + temp + ibh + dpg + ibt + vis,
#                 data = ozone[x_train_id, ], method = "BAS",
#                 prior = "hyper-g", alpha = 3, modelprior = uniform())
# 
#   ret2 = bas.lm(y ~ vh + wind + humidity + temp + ibh + dpg + ibt + vis,
#                 data = ozone[x_train_id, ], method = "BAS",
#                 prior = "g-prior", modelprior = uniform())
# 
#   ret_GP = GP.fit(ozone[x_train_id, 1:8], ozone[x_train_id, 'y'], a = 0, b = 0, jitter = 0)
# 
#   # summarize
#   data_all[[ith.rep]] = x_train_id
#   ret_all[[ith.rep]] = list(hyper.g = ret1, g = ret2, GP = ret_GP)
#   temp = lapply(ret_all[[ith.rep]], function(k) summarize(k, ozone, x_train_id))
# 
#   MSE_all[ith.rep, ] = do.call("c", lapply(temp, function(k) k$MSE))
#   RISK_all[ith.rep, ] = do.call("c", lapply(temp, function(k) k$RISK))
# 
# }
# 
# save(MSE_all, RISK_all, file = "MSE.RData")
# save(MSE_all, RISK_all, ret_all, data_all, file = "MSE_big.RData")
# t2 = proc.time()
# cat("time took: ", t2 - t1) # took less than 10min in Mac

### MSE of each model (useful in submitted paper -- rethink its role later on) ------
# set.seed(20161)
# N.rep = 100
# 
# x = dat$ozone$x.rescale
# y = dat$ozone$y
# n_total = length(y)
# p = ifelse(is.null(ncol(x)), 1, ncol(x))
# num.model = 2^p
# model.id = p2model.id(p)
# 
# 
# MSE_all_part2 = matrix(NA, nrow = N.rep, ncol = 12)
# 
# for (ith.rep in 1:N.rep)
# {
#   x_train_id = sample(1:n_total, round(0.5 * n_total))
#   y_train = y[x_train_id]
#   x_train = x[x_train_id, ]
#   
#   x_valid = x[-x_train_id, ]
#   y_valid = y[-x_train_id]
#   
#   y = as.vector(y)
#   n = length(y_train)
#   
#   ret_GP = GP.fit(x_train, y_train, a = 0, b = 0, jitter = 0)
#   
#   KL.all = matrix(NA, nrow = 6, ncol = 2^p) # all 2^p models
#   count = 0
#   for (g.type in 1:3){
#     for (g.sigma.method in 1:2){
#       count = count + 1
#       KL.all[count, ] = unlist(lapply(model.id$idx, function(k)
#         data2KL(X.j = cbind(1, x_train[, k]), y_train,
#                 g = Inf , GP = ret_GP, type = g.type,
#                 sigma.method = g.sigma.method, a = 0, b = 0)))
#     }
#   }
#   
#   summary.D_Bayes = sapply(1:6, function(k)
#     conv.summary.KL(KL.all[k, ], model.id, n = length(y_train), num.model))
#   colnames(summary.D_Bayes) = c("KL1", "KL1s", "KL2", "KL2s", "KL3", "KL3s")
#   
#   ret_MSE = rep(NA, 12)
#   for (ith.method in 1:6){
#     Z = rep(NA, num.model)
#     ret = summary.D_Bayes[, ith.method]
#     model.idx = as.integer(ret[1:num.model])
#     model.prob = ret[(num.model + 1):(2 * num.model)]
#     Z[model.idx] <- model.prob
#     
#     Z = exp(Z)
#     post.prob = Z
#     # post.prob = exp(post.prob - max(post.prob))
#     post.prob = post.prob / sum(post.prob)
#     flag_median = which(apply(model.id$m * post.prob, 2, sum) > 0.5)
#     flag_highest = model.id$idx[[which.max(post.prob)]]
#     
#     if (length(flag_highest) > 0){
#       lm_fit = lm(y_train ~ x_train[, flag_highest])
#     } else {
#       lm_fit = lm(y_train ~ 1)
#     }
#     
#     yhat = cbind(1, x_valid[, flag_highest]) %*% lm_fit$coefficients
#     ret_MSE[ith.method] = mean((y_valid - yhat)^2)
#     
#     if (length(flag_median) > 0){
#       lm_fit = lm(y_train ~ x_train[, flag_median])
#     } else {
#       lm_fit = lm(y_train ~ 1)
#     }
#     
#     yhat = cbind(1, x_valid[, flag_median]) %*% lm_fit$coefficients
#     ret_MSE[ith.method + 6] = mean((y_valid - yhat)^2)
#   }
#   MSE_all_part2[ith.rep, ] = ret_MSE
#   print(round(ret_MSE, 1)) # use "print" inside a loop if you want to print out outputs
# }
# save(MSE_all_part2, file = "MSE_all_part2.RData")

### exponential weighting using dpost and RT (both MSE and RISK) -----

load("MSE_big.RData")
N.rep = 100
x = dat$ozone$x.rescale
y = dat$ozone$y
n_total = length(y)
p = ifelse(is.null(ncol(x)), 1, ncol(x))
num.model = 2^p
model.id = p2model.id(p)

## calculate all predictions and fitted values 
yhat_all = as.list(1:N.rep)
for (ith.rep in 1:N.rep)
{
  x_train_id = data_all[[ith.rep]]
  y_train = y[x_train_id]
  x_train = x[x_train_id, ]
  
  x_valid = x[-x_train_id, ]
  y_valid = y[-x_train_id]
  
  y = as.vector(y)
  n = length(y_train)
  
  # calculate MSE and RISK for all sub-models per replication 
  yhat_all_temp = matrix(NA, num.model, length(y)) 
  for (i in 1:num.model){
    flag = model.id$id[[i]]
    if (length(flag) > 0){
      lm_fit = lm(y_train ~ x_train[, flag])
    } else {
      lm_fit = lm(y_train ~ 1)
    }
    
    yhat = cbind(1, x_valid[, flag]) %*% lm_fit$coefficients
    yhat_all_temp[i, x_train_id] = lm_fit$fitted.values
    yhat_all_temp[i, -x_train_id] = yhat
  }
  
  yhat_all[[ith.rep]] = yhat_all_temp
}


conv1 <- function(x_train_id, response = y){ # return RISK and MSE
  function(est){
    res = est - response 
    c(mean(res[x_train_id]^2), mean(res[-x_train_id]^2))
  }
}
logkernel2weight <- function(logkernel){
  rescale.kernel = logkernel - max(logkernel)
  weight = exp(rescale.kernel) / sum(exp(rescale.kernel))
  return(weight)
}

MSE_RISK = array(NA, c(4, 9, N.rep))
y = as.vector(y)
weights_all = list(1:N.rep)

extract_logweight <- function(ith.rep, x, y, data_all, ret_all){
  
  p = ifelse(is.null(ncol(x)), 1, ncol(x))
  num.model = 2^p
  model.id = p2model.id(p)
  
  x_train_id = data_all[[ith.rep]]
  y_train = y[x_train_id]; x_train = x[x_train_id, ]
  x_valid = x[-x_train_id, ]; y_valid = y[-x_train_id]
  n = length(y_train)
  
  ret_GP = ret_all[[ith.rep]]$GP
  
  KL.all = matrix(NA, nrow = 6, ncol = 2^p) # all 2^p models
  count = 0
  for (g.type in 1:3){
    for (g.sigma.method in 1:2){
      count = count + 1
      KL.all[count, ] = unlist(lapply(model.id$idx, function(k)
        data2KL(X.j = cbind(1, x_train[, k]), y_train,
                g = Inf , GP = ret_GP, type = g.type,
                sigma.method = g.sigma.method, a = 0, b = 0)))
    }
  }
  
  summary.D_Bayes = sapply(1:6, function(k)
    conv.summary.KL(KL.all[k, ], model.id, n = length(y_train), num.model))
  colnames(summary.D_Bayes) = c("KL1", "KL1s", "KL2", "KL2s", "KL3", "KL3s")
  
  log_weight = matrix(NA, num.model, 6)
  for (ith.method in 1:6){
    Z = rep(NA, num.model)
    ret = summary.D_Bayes[, ith.method]
    model.idx = as.integer(ret[1:num.model])
    model.prob = ret[(num.model + 1):(2 * num.model)]
    Z[model.idx] <- model.prob
    log_weight[,ith.method] = Z
  }
  
  return(log_weight)
}

justify <- function(bas.model, my.model){
  remove.intercept = bas.model[-1]
  if (length(remove.intercept) != length(my.model)){
    return(FALSE)
  } else {
    return(sum((remove.intercept - my.model)^2) == 0)
  }
}
for (ith.rep in 1:N.rep)
{
  ret_MSE = matrix(NA, 4, 9)
  x_train_id = data_all[[ith.rep]]
  
  Dpost = extract_logweight(ith.rep, x, y, data_all, ret_all)
  
  # use posterior mean in GP for sigma.square in weight.RT
  sigma.sq = ret_all[[ith.rep]]$GP$detail$sigma.sq
  beta = 4 * sigma.sq
  trace.Aj = as.numeric(lapply(model.id$id, length)) + 1
  RISK = as.numeric(apply(yhat_all[[ith.rep]], 1, conv1(x_train_id))[1,])
  Rn.unb = RISK + 2 * sigma.sq / length(x_train_id) * trace.Aj # up to a constant
  log.weight.RT =  -length(x_train_id) * Rn.unb - log(beta) # + log_prior_C
  
  # BMA
  BMA.ret = ret_all[[ith.rep]]$g
  my.idx = rep(NA, num.model)
  for (i in 1:num.model){
    my.idx[i] = which(as.numeric((lapply(BMA.ret$which, function(k) justify(k, model.id$idx[[i]])))) == 1)
  }
  BMA.weight = BMA.ret$logmarg[my.idx]
  
  BMAhyperg.ret = ret_all[[ith.rep]]$hyper.g
  my.idx = rep(NA, num.model)
  for (i in 1:num.model){
    my.idx[i] = which(as.numeric((lapply(BMAhyperg.ret$which, function(k) justify(k, model.id$idx[[i]])))) == 1)
  }
  BMAhyperg.weight = BMAhyperg.ret$logmarg[my.idx]
  
  log.weight = cbind(Dpost, log.weight.RT, BMA.weight, BMAhyperg.weight)
  colnames(log.weight) <- c("KL1", "KL1s", "KL2", "KL2s", "KL3", "KL3s", "RT", "BMA.g", "BMA.hyper.g")
  
  weights_all[[ith.rep]] = list(log.weight = log.weight, weight1 = 1, weight2 = 1)
  
  for (ith.method in 1:9){
    log_prior_C = rep(0, num.model) # equal weight prior
    weight = logkernel2weight(log.weight[,ith.method] + log_prior_C)
    weight_fit = colSums(yhat_all[[ith.rep]] * weight)
    ret_MSE[1:2, ith.method] = conv1(x_train_id)(weight_fit)
    weights_all[[ith.rep]]$weight1 = weight 
    
    M = num.model
    H_M = (1 - exp(-(M + 1))) / (1 - exp(-1))
    abs.p = as.numeric(lapply(model.id$id, length)) + 1
    log_prior_C = -(lchoose(M, abs.p) + abs.p + log(H_M))
    weight = logkernel2weight(log.weight[,ith.method] + log_prior_C)
    weight_fit = colSums(yhat_all[[ith.rep]] * weight)
    ret_MSE[3:4, ith.method] = conv1(x_train_id)(weight_fit)
    weights_all[[ith.rep]]$weight2 = weight 
  }
  
  
  MSE_RISK[,,ith.rep] = ret_MSE
}

colnames(MSE_RISK) = colnames(log.weight)
rownames(MSE_RISK) = c("RISK.uniform", "MSE.uniform", "RISK.binom", "MSE.binom") # different model prior

save(MSE_RISK, weights_all, file = "MSE_RISK.RData")

## summarize exponential weights ----
load("MSE_RISK.RData")
par(mfrow = c(1, 1))
boxplot(t(MSE_RISK["MSE.uniform", ,]))
boxplot(t(MSE_RISK["RISK.uniform", ,]))
boxplot(t(MSE_RISK["MSE.binom", ,]))
boxplot(t(MSE_RISK["RISK.binom", ,]))

boxplot(sqrt(t(MSE_RISK["MSE.uniform", c("KL1s", "KL2s", "RT", "BMA.g"),])))

# ESS
logweight2eff <- function(logweight){
  if (is.vector(logweight)){
    1/sum(logkernel2weight(logweight)^2)
  }
  if (is.matrix(logweight)){
    apply(logweight, 2, function(k) 1/sum(logkernel2weight(k)^2))
  }
  
}


#round 2 submission plot 
load("MSE_RISK_50.RData")
# plot ratio of RMSE vs RT 
# plot.data = sqrt(t(MSE_RISK["MSE.uniform", c("KL1s", "KL2s", "BMA.g"),])) / sqrt((MSE_RISK["MSE.uniform", c("RT"),]))
plot.data = sqrt(t(MSE_RISK["MSE.uniform", c("KL1s", "KL2s", "BMA.g"),])) - sqrt((MSE_RISK["MSE.uniform", c("RT"),]))
colnames(plot.data) = c("D1", "D2", "BMA")

png("RMSE.png")
boxplot(plot.data,cex.axis = 1.5, lwd = 1, las = 1)
abline(h = 0, lwd = 1.5, lty = 2)
dev.off()


a1 = sqrt((MSE_RISK[c("MSE.uniform","RISK.uniform"), c("KL1s", "KL2s", "RT", "BMA.g"),]))
ess1 = sapply(weights_all, function(k) logweight2eff(k$log.weight))
table.data = array(NA, c(3, 4, 100))
table.data[1:2, , ] = a1
table.data[3, , ] = ess1[c("KL1s", "KL2s", "RT", "BMA.g"), ]
table.data = apply(table.data, c(1,2), mean)
round(table.data, 2)
xtable(round(table.data, 3))


plot.data = t(ess1[c("KL1s", "KL2s", "RT", "BMA.g"), ])
png("ESS.png")
colnames(plot.data) = c("D1", "D2", "EW", "BMA")
boxplot(plot.data, cex.axis = 1.5, lwd = 1, las = 1)
dev.off()


### past -----

load("MSE.RData")
apply(sqrt(MSE_all), 2, mean)
sqrt(apply(sqrt(MSE_all), 2, var) / N.rep)

load("MSE_all_part2.RData")
apply(sqrt(MSE_all_part2), 2, mean)
sqrt(apply(sqrt(MSE_all_part2), 2, var) / N.rep)

# use KL1, KL2, g-prior, hyper-g: MSE_all_part2 col 1, col 3, MSE_all col 4, col 2
# reference model: MSE_all col 5
MSE_table = cbind(MSE_all_part2[, c(1, 3)], MSE_all[, c(4, 2)], MSE_all[, 5])
apply(sqrt(MSE_table), 2, mean)
sqrt(apply(sqrt(MSE_table), 2, var) / N.rep)
# > round(apply(sqrt(MSE_table), 2, mean), 2)
# [1] 4.61 4.61 4.62 4.63 4.09
# > max(sqrt(apply(sqrt(MSE_table), 2, var) / N.rep))
# [1] 0.01954885


## exponential weighting using dpost and RT (both MSE and prediction)

load("MSE_big.RData")
N.rep = 100
x = dat$ozone$x.rescale
y = dat$ozone$y
n_total = length(y)
p = ifelse(is.null(ncol(x)), 1, ncol(x))
num.model = 2^p
model.id = p2model.id(p)


MSE_all_part2 = matrix(NA, nrow = N.rep, ncol = 12)

for (ith.rep in 1:N.rep)
{
  x_train_id = sample(1:n_total, round(0.5 * n_total))
  y_train = y[x_train_id]
  x_train = x[x_train_id, ]
  
  x_valid = x[-x_train_id, ]
  y_valid = y[-x_train_id]
  
  y = as.vector(y)
  n = length(y_train)
  
  ret_GP = GP.fit(x_train, y_train, a = 0, b = 0, jitter = 0)
  
  KL.all = matrix(NA, nrow = 6, ncol = 2^p) # all 2^p models
  count = 0
  for (g.type in 1:3){
    for (g.sigma.method in 1:2){
      count = count + 1
      KL.all[count, ] = unlist(lapply(model.id$idx, function(k)
        data2KL(X.j = cbind(1, x_train[, k]), y_train,
                g = Inf , GP = ret_GP, type = g.type,
                sigma.method = g.sigma.method, a = 0, b = 0)))
    }
  }
  
  summary.D_Bayes = sapply(1:6, function(k)
    conv.summary.KL(KL.all[k, ], model.id, n = length(y_train), num.model))
  colnames(summary.D_Bayes) = c("KL1", "KL1s", "KL2", "KL2s", "KL3", "KL3s")
  
  ret_MSE = rep(NA, 12)
  for (ith.method in 1:6){
    Z = rep(NA, num.model)
    ret = summary.D_Bayes[, ith.method]
    model.idx = as.integer(ret[1:num.model])
    model.prob = ret[(num.model + 1):(2 * num.model)]
    Z[model.idx] <- model.prob
    
    Z = exp(Z)
    post.prob = Z
    # post.prob = exp(post.prob - max(post.prob))
    post.prob = post.prob / sum(post.prob)
    flag_median = which(apply(model.id$m * post.prob, 2, sum) > 0.5)
    flag_highest = model.id$idx[[which.max(post.prob)]]
    
    if (length(flag_highest) > 0){
      lm_fit = lm(y_train ~ x_train[, flag_highest])
    } else {
      lm_fit = lm(y_train ~ 1)
    }
    
    yhat = cbind(1, x_valid[, flag_highest]) %*% lm_fit$coefficients
    ret_MSE[ith.method] = mean((y_valid - yhat)^2)
    
    if (length(flag_median) > 0){
      lm_fit = lm(y_train ~ x_train[, flag_median])
    } else {
      lm_fit = lm(y_train ~ 1)
    }
    
    yhat = cbind(1, x_valid[, flag_median]) %*% lm_fit$coefficients
    ret_MSE[ith.method + 6] = mean((y_valid - yhat)^2)
  }
  MSE_all_part2[ith.rep, ] = ret_MSE
  print(round(ret_MSE, 1)) # use "print" inside a loop if you want to print out outputs
}

####################################################
##############Summary: uscrime #####################
####################################################
rm(list = ls())
load("ret.RData")

my.plot <- function(x, y, ...){
  upper = max(max(x), max(y))
  plot(x, y, xlim = c(0, upper), ylim = c(0, upper), ...)
  abline(c(0, 1))
}

## Effect of prior specification
## g = n gives ridicularly small D-prob. for small sample size 
## because everthing is in the exponential 
## For large sample size, g doesn't make any difference 
## # weird problem: Kl1s, KL2s, KL3s - all 0 - solved by rescaling 

x = dat$uscrime$x.rescale
y = dat$uscrime$y

p = ncol(x)
num.model = 2^p
method.name <- c("KL1", "KL1s", "KL2", "KL2s", "KL3", "KL3s", "g-prior", "hyper-g")

Z_all = list(1:2) # two priors where g = Inf vs g = n

ret = ret_uscrime_1; 
head(ret[num.model:(2 * num.model), ])
Z = matrix(NA, nrow = num.model, ncol = length(method.name))
colnames(Z) = method.name

for (i in method.name){
  model.id = as.integer(ret[1:num.model, i])
  model.prob = ret[(num.model + 1):(2 * num.model), i]
  Z[model.id, i] <- model.prob
}
Z[, 1:6] = exp(Z[, 1:6])
Z_all[[1]] = Z

ret = ret_uscrime_2; 
head(ret[num.model:(2 * num.model), ])
Z = matrix(NA, nrow = num.model, ncol = length(method.name))
colnames(Z) = method.name

for (i in method.name){
  model.id = as.integer(ret[1:num.model, i])
  model.prob = ret[(num.model + 1):(2 * num.model), i]
  Z[model.id, i] <- model.prob
}
Z[, 1:6] = exp(Z[, 1:6])
Z_all[[2]] = Z

# conditional D-bayes 
conditional_Z_all = lapply(Z_all, 
                           function(k) t(t(k[, 1:6]) / apply(k[, 1:6], 2, sum)))

# conditional D-bayes is much simillar than D-Bayes across type 1, 2, 3, priors 
my.plot(Z_all[[1]][, "KL1"], Z_all[[2]][, "KL1"])
my.plot(conditional_Z_all[[1]][, "KL1"], conditional_Z_all[[2]][, "KL1"])

my.plot(Z_all[[1]][, "KL1"], Z_all[[1]][, "KL2"])
my.plot(conditional_Z_all[[1]][, "KL1"], conditional_Z_all[[1]][, "KL2"])

## priors are influential - we use default priors g = Inf 
# # uscrime data - effect of priors 
# par(mfrow = c(1, 2))
# my.plot(Z_all[[1]][, "g-prior"], Z_all[[2]][, "hyper-g"], pch = 16, 
#         cex.axis = 2.5, lwd = 3, xlab = "g-prior", ylab = "hyper-g prior")
# 
# my.plot(Z_all[[1]][, "KL2"], Z_all[[2]][, "KL2"], cex.axis = 2.5, lwd = 3, pch = 16, 
#         xlab = "g = infinity", ylab = "g = n")

# uscrime data - conditional D-bayes | D-bayes  vs g-prior 

# idx = which(Z_all[[1]][, "g-prior"] > 1e-5)
idx = 1:nrow(Z_all[[1]])
my.plot(conditional_Z_all[[1]][idx, "KL2"], Z_all[[1]][idx, "g-prior"], pch = 16, 
        cex.axis = 2.5, lwd = 3, xlab = "Conditional D-bayes", ylab = "g prior")
my.plot(Z_all[[1]][idx, "KL2"], Z_all[[1]][idx, "g-prior"], pch = 16, 
        cex.axis = 2.5, lwd = 3, xlab = "D-bayes", ylab = "g prior")

# # uscrime data - conditional D-bayes | D-bayes  vs hyper-g-prior
# my.plot(conditional_Z_all[[1]][, "KL2"], Z_all[[1]][, "hyper-g"], pch = 16,
#         cex.axis = 2.5, lwd = 3, xlab = "Conditional D-bayes", ylab = "hyper-g prior")
# my.plot(Z_all[[1]][, "KL2"], Z_all[[1]][, "hyper-g"], pch = 16,
#         cex.axis = 2.5, lwd = 3, xlab = "D-bayes", ylab = "hyper-g prior")


# Not similar for uscrime data 
# # uscrime data - different D-bayes - different conditional D-bayes 
# my.plot(conditional_Z_all[[1]][, "KL1"], conditional_Z_all[[1]][, "KL2"], pch = 16, 
#         cex.axis = 2.5, lwd = 3, xlab = "Type 1", ylab = "Type 2")
# my.plot(conditional_Z_all[[1]][, "KL1"], conditional_Z_all[[1]][, "KL1s"], pch = 16, 
#         cex.axis = 2.5, lwd = 3, xlab = "Type 1", ylab = "Type 1 - simplified")


## selected model 
model.id = p2model.id(p)

post.prob = Z_all[[1]][, 'KL3s']
model.idx = which(rank(post.prob) > (2^p - 5))
model.id$idx[model.idx]

### nice scatterplot
### Ref: http://stats.stackexchange.com/questions/24416/how-to-do-a-pretty-scatter-plot-in-r
library(ggplot2) #make sure the newest is installed

df = data.frame(v1 = Z[,1], v2 = Z[,3])

# df <- data.frame(v1 = runif(1000), v2 = runif(1000))
# 
bin.plot<-qplot(data=df,
                x=v1,
                y=v2,
                z=v2)
bin.plot+stat_summary_hex(fun=function(z)length(z))


## summarize what we got 
g.ab = 1

n = length(y)
p = ncol(x)
m <- t(sapply(1:(2^p - 1),function(x){ as.integer(intToBits(x))[1:p]}))
idx = lapply(1:(2^p - 1), function(k) which(m[k, ] == 1))

ret = GP.fit(x, y, a = g.ab, b = g.ab, jitter = 0)
ret
mean((y - ret$detail$Y.hat)^2)

file.name = sprintf("dataApplication(ab = %.2f).RData", g.ab)
load(file.name)

KL0 = rep(NA, 4)
KL0.name = rep(NA, 4)
count = 0
for (g.type in 1:2){
  for (g.sigma.method in 1:2){
    count = count + 1
    KL0[count] = 
      data2KL(X.j = matrix(1, nrow = n, ncol = 1), y, g = Inf, GP = ret, type = g.type, 
              sigma.method = g.sigma.method, a = g.ab, b = g.ab)
    KL0.name[count] = sprintf('Type%d|sigma%d', g.type, g.sigma.method)
  }
}
names(KL0) = KL0.name

# add null model 
m = rbind(m, 0)
KL.all = cbind(KL.all, KL0)

abs.prob = exp(-n * KL.all)
rel.prob = abs.prob / rowSums(abs.prob)

# doesn't make sense so far 
inclusion = matrix(NA, nrow = 15, ncol = 4) # 4 methods 
for (j in 1:4){
  temp = m * c(rel.prob[j, ])
  inclusion[,j] = colSums(temp)
}
inclusion

previous = read.csv(file = "previous.csv", header = FALSE)

# include 10 models and calculate the relative prob and absolute prob 

indicator2idx <- function(indicator){
  sum(indicator * 2^(0:(p - 1)))
}
# m.recover = apply(m, 1, indicator2idx)
# all.equal(m.recover, c(1:(2^p - 1), 0)) # TRUE 

indicator.list = lapply(1:10, function(k) which(previous[, k + 1] > 0.5))
model.list = unique(apply(previous[, -1], 2, function(k) indicator2idx(k > 0.5)))
KL.list = KL.all[, model.list]
prob.rel.list = rel.prob[, model.list]
prob.list = exp(-n * KL.list)
round(prob.list, 3)

# conclude using top 5 models 
model.id = which(rank(abs.prob[4, ]) > (2^p - 5) )
abs.prob[4, model.id]
sapply(abs.prob[4, model.list], function(k) sum(abs.prob[4,] > k))

# new x adding nonlinear effects
target.id = idx[model.id][[2]]
conf <- function(target.id, new.x){
  g.type = 2; g.sigma.method = 2; 
  temp = data2KL(X.j = cbind(1, x[, target.id], new.x), y, g = Inf, GP = ret, 
                 type = g.type, sigma.method = g.sigma.method, a = g.ab, b = g.ab)
  exp(-n * temp)
}

# try goodness-of-fit
conf <- function(target.id, new.x){
  g.type = 2; g.sigma.method = 1; 
  temp = data2KL(X.j = cbind(1, x[, target.id], new.x), y, g = Inf, GP = ret, 
                 type = g.type, sigma.method = g.sigma.method, a = g.ab, b = g.ab)
  exp(-n * temp)
}
conf(c(1, 3, 4, 9,13,14), x[, 10]) # 0.5968
conf(c(1, 3, 4, 9,13,14), x[, 10]^8.5) # 0.6132
conf(c(1, 3, 4, 9,13,14), exp(x[, 10])) # 0.6085
conf(c(1, 3, 4, 9,13,14), exp(x[, 10] * 2)) # 0.6143
conf(c(1, 3, 4, 9,13,14), cos(1.5 * x[,10])) # 0.6197
conf(c(1, 3, 4, 9,13,14), cos(pi * 5/2 * x[,10])) # 0.6247
conf(c(1, 3, 4, 9,13,14), exp(cos(pi * 5/2 * x[,10])/2)) # 0.6284

X1 = x[,c(1,3,4,9,10,13,14)]
X2 = X1^2
X3 = cos(pi * 5/2 * X1)

dictionary = lapply(1:2^7, function(k) as.integer(intToBits(k))[1:7])
design = lapply(dictionary, function(k){
  X = X1; 
  if (length(which(k == 1)) > 0){
    X[,which(k == 1)] = X3[, which(k == 1)]
  }
  return(X)
})

test = lapply(design, function(k) data2KL(X.j = cbind(1, k), y, g = Inf, GP = ret, type = 2, 
                                          sigma.method = 2, a = g.ab, b = g.ab))
test = exp(-n * unlist(test)) # from 0.44 to 0.6006
summary(test)

## use b-spline for single variable and then select the best basis 
test = bs(x[,10], df = 15)
KL.all.test = matrix(NA, nrow = 4, ncol = 2^p - 1)
count = 0
for (g.type in 1:2){
  for (g.sigma.method in 1:2){
    count = count + 1
    KL.all.test[count, ] = unlist(lapply(idx, function(k) 
      data2KL(X.j = cbind(1, x[,c(1, 3, 4, 9,10, 13,14)], test[, k]), 
              y, g = n, GP = ret, type = g.type, 
              sigma.method = g.sigma.method, a = g.ab, b = g.ab)))
  }
}

abs.prob.test = exp(-n * KL.all.test)
min(abs.prob.test)
max(abs.prob.test)



