# Demonstration 

rm(list = ls())
source("sourceCode.R")

### Ozone data
library(earth)
data("ozone1")
y = ozone1[, 'O3']
x = as.matrix(ozone1[, -c(1, 10)]) 

# rescale x such that x \in [0, 1]^d 
rescale_x = x
for (i in 1:ncol(x)){
  rescale_x[,i] = (x[, i] - min(x[,i])) / (max(x[,i]) - min(x[,i]))
}
x = rescale_x

# fit the reference model 
g.ab = 0 
ret = GP.fit(x, y, a = g.ab, b = g.ab, jitter = 0)
# mean((y - ret$detail$Y.hat)^2) # residual sum of squares 

p = ncol(x)
m <- t(sapply(1:2^p,function(x){ as.integer(intToBits(x))[1:p]}))
idx = lapply(1:2^p, function(k) which(m[k, ] == 1))

# calculate Kullback-Leibler divergences for all sub-models 
# g = Inf means we use flat priors for coefficients
KL.all = matrix(NA, nrow = 2, ncol = 2^p) 
for (g.type in 1:2){
  KL.all[g.type, ] = unlist(lapply(idx, function(k)
    data2KL(X.j = cbind(1, x[, k]), y, g = Inf, GP = ret, type = g.type, a = g.ab, b = g.ab)))
}

n = length(y)
abs.prob = exp(-n * KL.all)
rel.prob = abs.prob / rowSums(abs.prob)

## chose the best model according to D-probabilities 
# use (type 1, type 2) D-probability 
idx[c(which.max(abs.prob[1, ]), which.max(abs.prob[2, ]))]

## use absolute D-probabilities to assess model fitting 
c(Type1=max(abs.prob[1,]), Type2=max(abs.prob[2,]))
# small D-probabilities suggest a lack of fit for all models 

## out-of-sample prediction error for selected models and the reference model 
set.seed(20161)
n = length(y)
x_train_id = sample(1:n, round(0.5 * n))
object = GP.fit(x[x_train_id, ], y[x_train_id], a = 0, b = 0, jitter = 0)
yhat_ref = predict(object, x_train = x[x_train_id, ],x_new = x[-x_train_id, ])
RMSE_ref = sqrt(mean((y[-x_train_id] - yhat_ref)^2))


x_select = x[, idx[[which.max(abs.prob[1, ])]]]
object = lm(y[x_train_id] ~ x_select[x_train_id, ])
yhat1 = cbind(1, x_select[-x_train_id, ]) %*% coef(object)
RMSE1 = sqrt(mean((y[-x_train_id] - yhat1)^2))

x_select = x[, idx[[which.max(abs.prob[2, ])]]]
object = lm(y[x_train_id] ~ x_select[x_train_id, ])
yhat2 = cbind(1, x_select[-x_train_id, ]) %*% coef(object)
RMSE2 = sqrt(mean((y[-x_train_id] - yhat2)^2))

# compare root MSE
c(Type1=RMSE1, Type2 = RMSE2, Reference=RMSE_ref)

save.image(file = "demo.RData")
