########### data generation #########
gen = function(p, n, d, coefs, seed = 2, rho = 0, intercept = 0, family = 'gaussian', ECcov = FALSE, tdf = NULL){
  y0 = NULL
	set.seed(seed)
	beta=c(coefs,rep(0,p-length(coefs)))

  cholmat.file = paste0('cholmat',p,'-rho',rho,'.txt')
	if(file.exists(cholmat.file)){
	  cholmat = read.table(file = cholmat.file)
	} else{
	  corrmat=rho^(abs(outer(1:p,1:p,'-')))
	  cholmat=chol(corrmat)
	  write.table(cholmat, file = cholmat.file)
	}
	
  cholmat = as.matrix(cholmat)
	x = matrix(rnorm(n*p),n,p)
	X = x%*%cholmat
	
	if(family == 'poisson'){
		lambda = exp(X%*%beta)
        y = rpois(n,lambda)
		}

	if(family == 'binomial'){
        mu = 1/(1+exp(-X%*%beta-intercept))
        y = rbinom(n,1,mu)
        
		}
	if(family == 'gaussian'){
	  if(is.null(tdf)){
	    eps = rnorm(n, 0, 1)
	  }else{
	    eps = rt(n, df = tdf)
	  }
		y0 = X %*% beta
		y = y0 + eps
		
		return(list(beta = beta, X = X, y = y, y0 = y0, eps = eps))
	    }
	
	return(list(beta = beta, X = X, y = y))
	}


### main function
f.main = function(b, data, nc, penalty = 'LASSO', family = 'gaussian'){
  
  ### data description
  n = dim(data$X)[1]
  p = dim(data$X)[2]
  
  ### compute the penalized path
  if(penalty == 'LASSO') fit.whole = glmnet(data$X, data$y, family = family)
  if(penalty != 'LASSO') {
    fit.whole = ncvreg(data$X, data$y, family = family, penalty = penalty)
    fit.whole$a0 = fit.whole$beta[1,]
    fit.whole$beta = fit.whole$beta[-1,]
    
  }
 
  dfs = apply(fit.whole$beta,2, f<-function(x){sum(x!=0)})
  lam.cand = which(dfs>0 & dfs < nc)
  model.seq.all = apply(fit.whole$beta[,lam.cand],2,f<-function(x){which(x!=0)})
  models = unique(model.seq.all)
  n.models = length(models)
  
  lam.list = fit.whole$lambda
  logliks = matrix(0,n.models, b)
  for(model.ind in 1:n.models){
    model = models[[model.ind]]
    for(reps in 1:b){
      set.seed(reps)
      cons.ind = sample(1:n)[1:nc]
      vali.ind = setdiff(1:n,cons.ind)
      
      fit.sub = glm(data$y[cons.ind] ~ data$X[cons.ind,model], family = family)
      
      logliks[model.ind, reps] = loglik(cbind(1,data$X[vali.ind,model]), data$y[vali.ind], fit.sub$coef, family = family)
    }
  }
  m.loglik = apply(logliks, 1, mean)
  m.ind = which.min(m.loglik) ###best model
  for(i in 1:length(model.seq.all)){
    if(setequal(model.seq.all[[i]],models[[m.ind]])){
      loc = i
      break
    }
  }
  lam.loc = lam.cand[loc]
  
  beta = fit.whole$beta[,lam.loc]
  a0 = fit.whole$a0[lam.loc]
  lambda = lam.list[lam.loc]
  
  return(list(model = models[[m.ind]], beta = beta, a0 = a0, lambda = lambda))
  
  
}


##log-likelihood
loglik <- function(X, y, beta, family) {
  link = as.vector(X %*% beta)
  n = length(y) 
  if (family == "gaussian") 
    return(n * log(mean((y - link)^2)))
  if (family == "poisson") 
    return(-2 * sum(exp(link) + 2 * y * link))
  if (family == "binomial") 
    return(2 * sum(log(1 + exp(link)) - y * link))
  
}


f.loss = function(y, yHat, family){
  
  if(family == 'gaussian'){
    def = (y - yHat)^2
    return(apply(def, 2, mean))
  }
  
  
  if(family == 'poisson'){
    if(length(devy[y == 0]) > 0){
      deveta = y * log(yHat) - yHat
      devy = y * log(y) - y
      devy[y == 0] = 0
      return(2*apply((devy-deveta),2,mean))
    }
    else return(0)
  }
  
  if(family == 'binomial'){
    return(1-apply(yHat, 2, f<-function(x){mean((x>0)==y)}))
    #return(apply(def, 2, mean))
  }
  
  
}



