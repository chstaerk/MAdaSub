library("MASS")
library("RcppNumerical")
#library("BMS")
library("BAS")
library("parallel")

library("mvnfast")
library("Rfast")


# MAdaSub algorithm (serial updating)

# Input: 
# data = (data$x,data$y) ("data" is list with data$x as the design matrix and data$y as the response vector)
# Iter: number of iterations
# priormean: initial proposal probabilities (prior mean of beta distribution). Can be the same (numeric) or different (vector) for all covariates
# L: parameters conrolling the adaptation of the algorithm. Can be the same (numeric, e.g. L=p) or different (vector) for all covariates
# family: "normal" for linear regression, "binomial" for logistic regression, "poisson" for poisson regression
# epsilon: "precision" parameter for truncation of proposal probabilities (e.g. epsilon=1/p)
# priorprob: prior probability of inclusion for the covariates (numeric)
# prior: specify the type of prior structure.
#        for family="normal" three different choices:
#        (1) set prior="independence" for use of independent prior on regression coefficients and Jeffreys prior on error variance
#        (2) set prior="gprior" for use of g-prior 
#        (3) set prior="EBIC" for use of prior corresponding to EBIC 
#        for family="binomial" or family="poisson" only prior choice (3) is available 
# g: constant for g-prior as well as independence prior (see also Remark 2.1 of paper)
# hyper: set hyper=TRUE for use of hyper beta-binomial model prior 
# a_prior: parameter a in beta-binomial model prior
# b_prior: parameter b in beta-binomial model prior
# const: constant "gamma" for prior corresponding to EBIC (gamma=0 for BIC)
# plot: if plot=TRUE, current trace plots of the posterior kernel are depicted (after every 10,000 iterations)
#       set plot=FALSE for computational speed-up
# S.start: optional initial model for MAdaSub sampler (otherwise randomly initialized)
# automatic_stop: set automatic_stop=TRUE for use of automatic stopping 
# precision: precision parameter for automatic stopping of the algorithm
# savings: only every "savings" iterations is saved (less memory usage for thinning)


# Output: 
# relfreq.hist: history of proposal probabilities (matrix of dimension p X floor(Iter/savings)) 
# relfreq.final: vector of final proposal probabilities (after the last iteration)
# importance.hist: history of empirical inclusion frequencies (matrix of dimension p X floor(Iter/savings)), 
# importance.final: final empirical inclusion frequencies 
# V.size, S.size: sizes of proposed models V and sampled models S in each iteration (vectors of dimension Iter)
# values: values of the log-kernel of the posterior of sampled models along the iterations (vector of dimension Iter)
# best.S: "Best" model (with the largest posterior probability) that has been found during MAdaSub  
# best.models: List of sampled models (of length Iter)
# acc.prob: Vector (of length Iter) of current acceptance rates along the iterations  
# acc: Absolute number of accepted Metropolis-Hastings moves
MAdaSub <-function (data,Iter,priormean,L,const=0,savings=1,family="normal",epsilon=1e-06,priorprob,prior,g=NULL,hyper=FALSE,a_prior=NULL,b_prior=NULL,plot=TRUE,S.start=NULL, automatic.stop =FALSE, precision = 0.001, burnin=0, s_max=NULL) { 
  
  p=ncol(data$x)
  n=nrow(data$x)
  
  if (is.null(s_max) & family =="normal" & prior=="EBIC") s_max = n-2
  if (is.null(s_max) & family =="normal" & prior=="gprior") s_max = n-1
  if (is.null(s_max) & family !="normal" & prior=="EBIC") s_max = n-1
  if (is.null(s_max)) s_max = p
  
  relfreq=numeric(p)
  relfreq[1:p] = priormean
  
  CS=matrix(NA,p,floor(Iter/savings))
  
  CS.cur = numeric(p) #numeric(p)
  CS.burnin = numeric(p) 
  
  V.size<- rep(NA,Iter) #numeric(Iter)
  S.size<- rep(NA,Iter) #numeric(Iter)
  
  values <- rep(NA,Iter)
  
  best.models = vector("list", Iter)
  
  if (is.null(S.start)) { 
    b = rbinom(p,1,relfreq)
    S = which(b==1)
    #S = integer(0)
  } else {
    S = S.start
  }
  
  acc.prob = rep(NA,Iter)  #numeric(Iter)
  acc = 0
  acc.after.burnin = 0
  
  for (t in 1:Iter) {
    
    relfreq[relfreq>1-epsilon] = 1-epsilon
    relfreq[relfreq<epsilon] = epsilon
    
    b = rbinom(p,1,relfreq)
    V = which(b==1)
    
    if (length(V)<=s_max)  { # reject V if length(V)>s_max
      if (prior=="EBIC"){
        alpha = (-1/2)* EBIC_glm(data=data,indices=V,const=const,family=family) +
          log_prob(indices=S,p=p,relfreq=relfreq) - 
          (-1/2)* EBIC_glm(data=data,indices=S,const=const,family=family) -
          log_prob(indices=V,p=p,relfreq=relfreq) 
      }   else {
        alpha =  marginal_log_like(data=data,indices=V,family=family,prior=prior,g=g) + 
          log_modelprior(indices=V,p=p,priorprob=priorprob,hyper=hyper,a_prior=a_prior,b_prior=b_prior) +
          log_prob(indices=S,p=p,relfreq=relfreq) - 
          marginal_log_like(data=data,indices=S,family=family,prior=prior,g=g)  - 
          log_modelprior(indices=S,p=p,priorprob=priorprob,hyper=hyper,a_prior=a_prior,b_prior=b_prior) -
          log_prob(indices=V,p=p,relfreq=relfreq)
      }
      u=log(runif(1))
      
      if (u<=alpha) { 
        S = V 
        acc = acc + 1
        if (t > burnin) {
          acc.after.burnin = acc.after.burnin + 1
        }
      }
    }
    
    
    if (prior=="EBIC"){
      values[t] = (-1/2) * EBIC_glm(data=data,indices=S,const=const,family=family) } else {
        values[t] = log_modelposterior(data=data,indices=S,p=p,priorprob=priorprob,family=family,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior) 
      }
    
    CS.cur[S] = CS.cur[S] +1
    V.size[t]=length(V)
    S.size[t]=length(S)
    
    relfreq = (L*priormean + CS.cur)/(L + t)
    
    if (t %% savings == 0) {
      CS[,t/savings] = CS.cur
    }
    
    if (t==burnin) CS.burnin <- CS.cur
    
    if (values[t]==max(values[1:t])) best.S = S
    
    best.models[[t]] = S
    
    acc.prob[t] = acc/t
    
    if (automatic.stop && all(abs(relfreq - CS.cur/t) <= precision )) break
    
    if (plot && (t %% 10000 == 0)) { 
      par(mfrow=c(1,1))
      plot(values,pch=20,main="MAdaSub" ,xlab="Iterations",ylab="Posterior")
      #abline(h=EBIC_glm(data,which(relfreq>0.5),const,family=family),col="blue")
      #abline(h=EBIC_glm(data,which(relfreq>0.9),const,family=family),col="red")
    }
  }
  
  best.S = best.models[[which.max(values)]]
  
  helpmat = matrix(rep((1:floor(Iter/savings))*savings,p) , p, floor(Iter/savings),byrow=TRUE)
  
  
  importance.hist = CS/helpmat 
  
  importance.final = (CS.cur-CS.burnin)/(Iter-burnin)
  #importance.final = CS.cur/Iter
  
  acc.prob.after.burnin <- acc.after.burnin / (Iter-burnin)
  
  relfreq.hist=(L*priormean + CS)/(L+helpmat) 
  relfreq.final = relfreq
  #relfreq.final=relfreq.hist[,floor(Iter/savings)]
  
  
  par(mfrow=c(1,1))
  if (plot) plot(values,pch=20,main="MAdaSub" ,xlab="Iterations",ylab="Posterior")
  #abline(h=EBIC_glm(data,which(relfreq>0.5),const,family=family),col="blue")
  #abline(h=EBIC_glm(data,which(relfreq>0.9),const,family=family),col="red")
  
  return(list(relfreq.hist=relfreq.hist, 
              relfreq.final=relfreq.final,
              importance.hist = importance.hist, 
              importance.final=importance.final,
              S.size=S.size,
              V.size=V.size,
              values=values,
              best.S=best.S,
              best.models=best.models, 
              acc.prob= acc.prob, 
              acc=acc, 
              CS=CS, 
              CS.cur=CS.cur, 
              S=S,
              acc.prob.after.burnin=acc.prob.after.burnin))
}

# MAdaSub algorithm with parallel updating

# Input (see also serial version above):
# nb.k: number of chains 
# M: number of iterations (T) per round for each chain 
# repeats: number of repeats 
# frac.par.update: fraction of chains with parallel updating (if frac.par.update=1, then all chains communicate)  
# burnin_rounds: number of burn-in rounds
# save_workers: number of chains (workers) for which the whole history of sampled models and acceptance rates is saved
#               (default: save_workers=1 to save memory)              

# Output (see also serial version above):
# parameter.list: List (of lenght nb.k) of final parameters for each MAdaSub chains 
#                 List of lists with components: priormean.cur as final proposal probabilities (of length p)
#                                                L.cur as final "precision" parameter
#                                                S.cur as final sampled model 
#                                                CS.cur as final absolute frequency each variable was sampled (of length p)
# importance.individual: 3-dimensional array with empirical inclusion frequencies of dimensions (repeats+1) X nb.k X p 
#                        First component: number of round (repeats+1 corresponds to final inclusion frequencies excluding burn-in rounds)
#                        Second component: number of chain 
#                        Third component: covariate index 
# acc: acceptance rates  
# models_parallel: sampled models for chains with parallel updating 
# models_serial: sampled models for chains with serial updating 
MAdaSub_parallel <- function (data,nb.k=2,M=10000,repeats=10,priormean,L,const=0,savings=1,family="normal",epsilon=1e-06,priorprob,prior,g=NULL,hyper=FALSE,a_prior=NULL,b_prior=NULL,
                              frac.par.update=1, burnin_rounds = 1, save_workers = 0, numCores = 1) { 
  p=ncol(data$x)
  
  parameter.list = vector("list", length = nb.k)
  for (l in 1:nb.k) {
    parameter.list[[l]]$CS.cur = numeric(p)
    parameter.list[[l]]$S.cur = NULL
  }
  
  if (length(L)==1) {
    for (l in 1:nb.k) {
      parameter.list[[l]]$L.cur = L
    }
  }
  
  if (length(L)==nb.k) {
    for (l in 1:nb.k) {
      parameter.list[[l]]$L.cur = L[[l]]
    }
  }
  
  if (length(priormean)==1) {
    for (l in 1:nb.k) {
      parameter.list[[l]]$priormean.cur = priormean
    }
  }
  
  if (length(priormean)==nb.k) {
    for (l in 1:nb.k) {
      parameter.list[[l]]$priormean.cur = priormean[[l]]
    }
  }
  
  
  # empirical inclusion frequencies 
  #(repeats+1 corresponds to overall incl. freq. excluding the first "burnin_rounds")
  importance.individual = array(NA, dim=c(repeats+1, nb.k, p))  
  
  # acceptance rates
  acc = array(NA, dim = c(repeats, nb.k)) 
  
  # sampled models
  models_parallel = rep(list(list()), save_workers)
  models_serial = rep(list(list()), save_workers)
  
  
  for (r in 1:repeats){
    
    results.cur = mclapply(parameter.list, f_parallel, mc.cores = numCores, 
                           data=data ,M=M ,const=const ,savings=savings ,family=family ,epsilon=epsilon ,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,
                           mc.set.seed = TRUE,  mc.preschedule = TRUE )
    skip.streams(3) # continue with random seeds
    
    #results[[r]] <- results.cur
    
    if(frac.par.update>0){  # chains with joint updating
      CS.cur.total = 0
      for (l in 1:(frac.par.update*nb.k)) {
        CS.cur.total = CS.cur.total + results.cur[[l]]$listout$CS.cur 
      }
      
      
      for (l in 1:(frac.par.update*nb.k)) {
        importance.individual[r,l,] = results.cur[[l]]$listout$CS.cur / M
        parameter.list[[l]]$S.cur = results.cur[[l]]$listout$S.cur
        
        # global updates
        parameter.list[[l]]$priormean.cur = (parameter.list[[l]]$L.cur*parameter.list[[l]]$priormean.cur + CS.cur.total) / ( parameter.list[[l]]$L.cur + M*nb.k*frac.par.update)
        parameter.list[[l]]$L.cur = parameter.list[[l]]$L.cur + M*nb.k*frac.par.update
        
        # save acceptance and sampled models
        acc[r, l] <- results.cur[[l]]$listout$acc
        if (l <= save_workers) {
          models_parallel[[r]][[l]] <- results.cur[[l]]$listout$best.models
        }
      }
    }
    
    if(frac.par.update<1){ # chains without joint updating
      for (l in (frac.par.update*nb.k+1):nb.k) {
        importance.individual[r,l,] = results.cur[[l]]$listout$CS.cur / M
        parameter.list[[l]]$S.cur = results.cur[[l]]$listout$S.cur
        
        # local updates
        parameter.list[[l]]$priormean.cur = (parameter.list[[l]]$L.cur*parameter.list[[l]]$priormean.cur + results.cur[[l]]$listout$CS.cur ) / ( parameter.list[[l]]$L.cur + M)
        parameter.list[[l]]$L.cur = parameter.list[[l]]$L.cur + M
        
        # save acceptance and sampled models
        acc[r, l] <- results.cur[[l]]$listout$acc
        if (l <= save_workers) {
          models_serial[[r]][[l]] <- results.cur[[l]]$listout$best.models
        }
      }
    }
  }
  importance.individual[repeats+1,,] = apply(importance.individual[-c(1:burnin_rounds,repeats+1),,], c(2,3), mean ) 
  return(list(parameter.list = parameter.list, importance.individual = importance.individual, acc = acc, models_parallel = models_parallel, models_serial = models_serial))
}

# Function to be parallelized (apply MAdaSub on different chains)
f_parallel <- function(parameter.list,data,M,const,savings,family,epsilon,priorprob,prior,g,hyper,a_prior,b_prior){
  output.cur = MAdaSub(data=data,Iter=M,priormean=parameter.list$priormean.cur ,L=parameter.list$L.cur,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,S.start=  parameter.list$S.cur,
                       plot=FALSE)
  listout = list()
  listout$S.cur = output.cur$S
  listout$CS.cur = output.cur$CS.cur
  listout$best.models = output.cur$best.models
  listout$acc = output.cur$acc
  
  return(list(listout=listout)) 
}


# Takes vector Iter = (Iter1,Iter2,...) and successively runs MAdaSub with Iter1, Iter2,... iterations
# using the updated priormean and L parameters. 
MAdaSub_wrapper <-function (data,Iter,priormean,L,const=0,savings=1,family="normal",epsilon=1e-06,priorprob,prior,g=NULL,hyper=FALSE,a_prior=NULL,b_prior=NULL, plot=TRUE) { 
  
  L.cur = L
  priormean.cur = priormean
  
  #if (save_all_output) output = list()
  
  for (i in 1:length(Iter)){
    
    start.time.cur <- Sys.time()
    output.cur = MAdaSub(data=data,Iter=Iter[i],priormean=priormean.cur,L=L.cur,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior, plot = plot)
    end.time.cur <- Sys.time()
    
    priormean.cur = output.cur$relfreq.final
    L.cur = L.cur + Iter[i]
    
    if (i==1) {
      print(paste("Iterations: ", 1, "-", Iter[i]))
    } else {
      print(paste("Iterations: ", sum(Iter[1:(i-1)])+1, "-", sum(Iter[1:i])))
    }
    print(paste("Acceptance rate: ", output.cur$acc/Iter[i]))
    print(end.time.cur - start.time.cur)
    #if (save_all_output)  output[[i]] = output.cur
  }
  #if (save_all_output) return(output) else return(output.cur)
  return(output.cur)
}



EBIC_glm <- function(data,indices,const,family="normal", maxiter = 20, precision =1e-02, start = NULL){
  n=nrow(data$x)
  p=ncol(data$x)
  x.cur=data$x[,indices]
  x.cur=cbind(c(rep(1,n)),x.cur)                 #include intercept!
  if (family=="normal") glm.out = .lm.fit(x.cur, data$y) # start=start) #,tol=2e-07) # intercept=TRUE ?
  if (family=="binomial") glm.out = fastLR(x.cur, data$y)
  if (family=="poisson") glm.out = glm.fit(x.cur, data$y, family=poisson(),intercept=FALSE,control = list(maxit = 20, epsilon=precision))# ,start=start) # intercept=TRUE ?
  if (family=="normal") {
    deviance = n*(1+log(2*pi)+ log(sum(glm.out$residuals^2) /n))
    EBIC = deviance + log(n)*(length(indices)+2) + 2*(length(indices)+2)*const*log(p)
  } else{
    deviance = -2 * glm.out$loglikelihood
    EBIC = deviance + log(n)*(length(indices)+1) + 2*(length(indices)+1)*const*log(p)}
  beta =glm.out$beta
  return(EBIC=EBIC)
}


# g-prior with g=n with Jeffrey prior on sigma^2 (unit information prior)
# or prior="independence" with constant g.
marginal_log_like <- function(data,indices,family="normal",prior="gprior",g){
  #n = nrow(data$x)
  #s = length(indices)
  #y = data$y
  #x.cur = data$x[,indices]
  #ybar = mean(data$y)
  #if (s>0){ m = -s/2 * log(g) - 1/2*det(t(x.cur) %*% x.cur + 1/g* diag(s)) - (n-1)/2 * log( sum((y-ybar)^2) -
  #    t(y)%*% x.cur %*% solve( t(x.cur) %*% x.cur + 1/g*diag(s)) %*% t(x.cur) %*% y  ) }else{
  #      m= - (n-1)/2 * log( sum((y-ybar)^2) - t(y)%*% y)   }
  #return(m)
  if (prior=="gprior") {
    s <- length(indices)
    fitted_lm <- lmfit(cbind(1, data$x[,indices]), data$y)
    R_squared <- 1 - sum(fitted_lm$residuals^2) / sum((data$y-mean(data$y))^2)
    marginal <- (n-s-1)/2 * log(1+g) - (n-1)/2 * log((1+g*(1-R_squared)))
    return(marginal)
  }
  #if (prior=="gprior") {
  #  data.cur = data
  #  data.cur$x = data$x[,indices]
  #  dataf = as.data.frame(data.cur)
  #  out = zlm(y~.,data=dataf)
  #  return(out$marg.lik)
  #}
  if (prior=="independence"){
    s = length(indices)
    n = length(data$y)
    if (s>0) {
      Vs_inv = 1/g * diag(s)
      data_x = data$x[,indices]
      Mridge = t(data_x) %*% data_x + Vs_inv
      Mridge_inv = ginv(Mridge)
      y_stan = data$y - mean(data$y)
      term_problematic = t(y_stan) %*%  y_stan - t(data$y) %*% data_x %*% Mridge_inv %*% t(data_x) %*% data$y
      if (term_problematic<0) {
        warning("Problems with marginal likelihood.")
        term_problematic=1
      }
      out = -1/2 * log(det(Mridge)) - 1/2 * s * log(g) - (n-1)/2 * log(term_problematic)
    } else {
      y_stan = data$y - mean(data$y)
      out = - (n-1)/2 * log( t(y_stan) %*%  y_stan )  
    }
    return(out)
  }
}


log_prob <- function(indices,p,relfreq){
  return( sum(log(relfreq[indices])) + sum(log(1-relfreq[(1:p)[-indices]])) )
}

log_modelprior <- function(indices,p,priorprob,hyper=FALSE,a_prior=NULL,b_prior=NULL){
  s = length(indices)
  if (hyper) { return(log(beta(a_prior + s , b_prior + p - s)) - log(beta(a_prior, b_prior)) )   } else {
    return(s*log(priorprob) + (p-s)*log(1-priorprob))
  }
}

log_modelposterior <- function(data,indices,p,priorprob,family="normal",prior,g,hyper=hyper,a_prior=a_prior,b_prior=b_prior){
  return( log_modelprior(indices=indices,p=p,priorprob=priorprob,hyper=hyper,a_prior=a_prior,b_prior=b_prior) + marginal_log_like(data=data,indices=indices,family=family,prior,g))
}

proposal_Griffin <- function(current, prop, A, D){
  return ( sum(log(A[current==0 & prop ==1])) +  sum(log(1-A[current==0 & prop ==0])) +
             sum(log(D[current==1 & prop ==0])) +  sum(log(1-D[current==1 & prop ==1]))  )
}




simdata <- function (n,p,beta,corr=0,family="normal",sigma.normal=1,corr.type="global",blocks=10) { # "blocks" = no. of blocks
  
  x=matrix(numeric(n*p),n,p) #initialize design matrix
  y=numeric(n) #initialize response variable
  
  mu=rep(0,p)
  
  if (corr.type=="global"){
    Sigma=diag(1-corr,p)+matrix(corr,p,p)
  }
  
  if (corr.type=="block"){
    Sigma=matrix(0,p,p)
    for (k in 1:p){
      for (m in 1:p){
        if ((k-m) %% blocks == 0 ) Sigma[k,m]= corr 
      }
    }
    Sigma=diag(1-corr,p)+Sigma
  }
  
  if (corr.type=="toeplitz"){
    help=numeric(p)
    for (k in 1:p) help[k]=corr^(k-1)
    Sigma=toeplitz(help)
  }
  
  if (p<=200) x = mvrnorm(n , mu, Sigma) else x = rmvn(n, mu, Sigma)
  #x = mvrnorm(n , mu, Sigma)
  
  linpred=x%*%beta
  
  if (family=="normal"){
    y=rnorm(n,linpred,sigma.normal)
  }
  
  if (family=="binomial"){
    #prob=exp(linpred)/(1+exp(linpred))
    prob = 1 / (1 + exp(-linpred))
    y=rbinom(n, size= 1, prob=prob)
  }
  
  if (family=="poisson"){
    lambda=exp(linpred)
    y=rpois(n, lambda=lambda)
  }
  
  return(list(x=x,y=y))
}

## helper for random seeds in parallel version
skip.streams <- function(n) {
  x <- .Random.seed
  for (i in seq_len(n))
    x <- nextRNGStream(x)
  assign('.Random.seed', x, pos=.GlobalEnv)
}


