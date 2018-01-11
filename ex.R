library('MASS')
library('glmnet')
library('ncvreg')

source('funs.R') ##load source file
family = 'gaussian'
nc.rate = 0.5 ##nc rate


n = 500 ##sample size
coefs = seq(from = 1,to = 0.4,by = -0.1)
 
p = 1000 ##no. of variable, can increase to 10000

b = 50 ##no. of random splits

oracle.ind = which(coefs!=0)
d = length(oracle.ind)
rho = 0
set.seed(0)
data = gen(p, n, d, coefs = coefs, seed = randSeed, rho = rho, family = family)
data.test = gen(p, n, d, coefs = coefs, seed = randSeed+1000, rho = rho, family = family) 
pen = 'LASSO'

nc = ceiling(n^nc.rate)

CCV.fit = f.main(b, data, nc, penalty = pen, family = family)
CCV.ind = which(CCV.fit$beta != 0)
CCV.FP = length(which(CCV.fit$beta[-oracle.ind] != 0))
CCV.FN = length(which(CCV.fit$beta[oracle.ind] == 0))

re.fit = glm(data$y ~ data$X[, CCV.ind],family=family)
re.yHat = t(t(as.matrix(data.test$X[, CCV.ind]) %*% re.fit$coefficients[-1]) + re.fit$coefficients[1])
CCV.PE = f.loss(data.test$y, re.yHat, family = family)

CCV.result = c(CCV.FP, CCV.FN, CCV.PE) ##a vector of FP, FN and prediction error on the test set



