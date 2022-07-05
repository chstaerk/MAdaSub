
library("BAS")
library("coda")

library("glmnet")

#RNGkind(kind = "default", normal.kind = "default", sample.kind = "default")
RNGkind("L'Ecuyer-CMRG")

#setwd("C:\\Users\\Staerk\\scieboBonn\\MAdaSub PAPER\\R files")
source("MAdaSub_LOAD.R")

# MC3 with option for swapping (set swap=TRUE)
MC3 <-function (Iter,data,priormean,L,const=0,savings=1,family="normal",epsilon=1e-06,priorprob,prior,g,hyper=FALSE,a_prior=NULL,b_prior=NULL,plot=TRUE, swap=FALSE, burnin=0, S.start=integer(0)) { 
  
  p=ncol(data$x)
  n=nrow(data$x)
  
  CV=matrix(numeric(Iter/savings*p),p,floor(Iter/savings))
  CS=matrix(numeric(Iter/savings*p),p,floor(Iter/savings))
  
  CV.cur = numeric(p)
  CS.cur = numeric(p)
  CS.burnin = numeric(p) 
  
  V.size<-numeric(Iter)
  S.size<-numeric(Iter)
  
  values <- rep(NA,Iter)
  
  best.models = vector("list", Iter)
  
  if (is.null(S.start)) { 
    b = rbinom(p,1,priormean)
    S = which(b==1)
    #S = integer(0)
  } else {
    S = S.start
  }
  
  if (swap) {
    sw <- rbinom(Iter, 1, 0.5) # Swapping move with probability 0.5, else flip (add/delete)
  }
  
  acc.prob = numeric(Iter)
  acc = 0
  acc.after.burnin = 0
  
  for (t in 1:Iter) {
    
    if (swap & sw[t]==1 & length(S)>0 & length(S)<p) { # Swapping 
      j1 = sample(S, 1)
      j2 = sample(setdiff(1:p, S), 1)
      V = union(setdiff(S, j1), j2) # Exclude j1 and include j2 (Swap)
    } else { # Flipping (add/delete)
      j = sample(1:p,1)
      V = S
      if (is.element(j,S)) {
        V = setdiff(V, j)
      } else {
        V = union(V, j)
      }
    }
    
    #if (length(V)>U_C){ V=sample(V,U_C) } 
    
    if (prior=="EBIC"){
      alpha = (-1/2)* EBIC_glm(data=data,indices=V,const=const,family=family)  - 
        (-1/2)* EBIC_glm(data=data,indices=S,const=const,family=family) 
    } else {
      alpha =  marginal_log_like(data=data,indices=V,family=family,prior=prior,g=g) + 
        log_modelprior(indices=V,p=p,priorprob=priorprob) - 
        marginal_log_like(data=data,indices=S,family=family,prior=prior,g=g)  - 
        log_modelprior(indices=S,p=p,priorprob=priorprob) 
    }
    u=log(runif(1))
    
    if (u<=alpha) { 
      S = V 
      acc = acc + 1
      if (t > burnin) {
        acc.after.burnin = acc.after.burnin + 1
      }
    }
    
    if (prior=="EBIC"){
      values[t] = (-1/2) * EBIC_glm(data=data,indices=S,const=const,family=family) 
    } else {
      values[t] = log_modelposterior(data=data,indices=S,p=p,priorprob=priorprob,family=family,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior) 
    }
    
    CV.cur[V] = CV.cur[V] +1
    CS.cur[S] = CS.cur[S] +1
    V.size[t]=length(V)
    S.size[t]=length(S)
    
    #relfreq = (L*priormean + CS.cur)/(L + t)
    
    if (t %% savings == 0) {
      CV[,t/savings] = CV.cur
      CS[,t/savings] = CS.cur
    }
    
    if (t==burnin) CS.burnin <- CS.cur
    
    if (values[t]==min(values[1:t])) best.S = S
    
    best.models[[t]] = S
    
    acc.prob[t] = acc/t
    
    if ((t %% 10000 == 0) && plot){ 
      #print(t)
      par(mfrow=c(1,1))
      plot(values,pch=20,main="MC3" ,xlab="Iter",ylab="value")
    }
  }
  
  best.S = best.models[[which.min(values)]]
  
  helpmat = matrix(rep((1:floor(Iter/savings))*savings,p) , p, floor(Iter/savings),byrow=TRUE)
  
  importance.hist = CS/helpmat 
  
  importance.final = (CS.cur-CS.burnin)/(Iter-burnin)
  #importance.final = CS.cur/Iter
  
  acc.prob.after.burnin <- acc.after.burnin / (Iter-burnin)
  
  par(mfrow=c(1,1))
  if (plot){
    plot(values,pch=20,main="MC3" ,xlab="Iter",ylab="value")
  }
  
  return(list(importance.final=importance.final,importance.hist=importance.hist,S.size=S.size,V.size=V.size,values=values,best.S=best.S,best.models=best.models, acc.prob= acc.prob, acc=acc, CS=CS,
              acc.prob.after.burnin = acc.prob.after.burnin))
}


##################################################################
# Simulate data as in Griffin et al. (2021)
##################################################################

family = "normal"
n = 500 
p = 500
SNR = 0.5   # SNR = 0.5: median r hat = 15.9 for 50000 iterations (n=500, p=500)
            # SNR = 2: median r hat = 7.9 for 200000 iterations (with 100000 burnin) (n=500, p=500)
sigma.normal = 1 

beta1=numeric(p)
S0=1:10
beta1[1:10]=c(2 , -3 , 2, 2, -3, 3, -2, 3, -2, 3) * (sigma.normal^2*log(p)/n)^(1/2) * SNR
S0=10
corr=0.6  
blocks=10
corr.type="toeplitz"
set.seed(21)
data=simdata(n,p,beta=beta1,corr= corr,family=family,sigma.normal=sigma.normal,blocks=blocks, corr.type=corr.type)

#### 

#lassofit <- cv.glmnet(data$x, data$y)
#lasso_estimate <- coef(lassofit)

##################################################################
# Begin MAdaSub algo
##################################################################

savings = 500 

q=10 

epsilon = 1/p ### ATTENTION epsilon = 0.1/p
priormean = q/p 
L = p 
priorprob= q/p 
g=9 

M <- 1000 
burnin_rounds <- 25 # 25
repeats <- 50 # 50 
Iterations <- M * repeats  
burnin <- M * burnin_rounds 
nb.k <- 5 
numCores <- 5 # nb.k

prior="independence" # g-prior (original form with Jeffreys prior on variance)
hyper = FALSE
a_prior = NULL
b_prior = NULL

swap = TRUE

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)

N_sim <- 200 # 200

Times_MAdaSub <- rep(NA,N_sim)
Times_MAdaSub_parallel <- rep(NA,N_sim)
Times_MAdaSub_parallel_without <- rep(NA,N_sim)
Times_MAdaSub_parallel_only_rj <- rep(NA,N_sim)
Times_MC3 <- rep(NA,N_sim)

Estimates_MAdaSub <- matrix(NA, N_sim, p)
Estimates_MAdaSub_parallel <- matrix(NA, N_sim, p)
Estimates_MAdaSub_parallel_without <- matrix(NA, N_sim, p)
Estimates_MAdaSub_parallel_only_rj <- matrix(NA, N_sim, p)
Estimates_MC3 <- matrix(NA, N_sim, p)

Accrate_MAdaSub <- rep(NA,N_sim)
Accrate_MAdaSub_parallel <- rep(NA,N_sim)
Accrate_MAdaSub_parallel_without <- rep(NA,N_sim)
Accrate_MAdaSub_parallel_only_rj <- rep(NA,N_sim)
Accrate_MC3 <- rep(NA,N_sim)

for (i in 1:N_sim){
  print(paste("Run ", i))
  set.seed(i)
  
  ###############################################
  ## with random initialization
  set.seed(i)
  ###
  priormean_parallel <- list()
  for (k in 1:nb.k) {
    priormean_parallel[[k]] <- rep(runif(1,2,10)/p, p)
  }
  L_parallel <- list()
  for (k in 1:nb.k) {
    L_parallel[[k]] <- runif(1, p/2, 2*p)
  }
  #L = p
  ###
  
  ## MAdaSub with parallel updating (frac.par.update=1)
  start.time <- Sys.time()
  output_parallel = MAdaSub_parallel(data=data,nb.k=nb.k,M=M,repeats=repeats,priormean=priormean_parallel,L=L_parallel,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=FALSE,a_prior=NULL,b_prior=NULL,
                             frac.par.update=1, burnin_rounds=burnin_rounds, numCores = numCores) 
  end.time <- Sys.time()
  Times_MAdaSub_parallel[i] <- difftime(end.time, start.time, unit = "sec")
  Estimates_MAdaSub_parallel[i,] <-   apply(output_parallel$importance.individual[repeats+1,,], 2, mean)
  Accrate_MAdaSub_parallel[i] <- mean(output_parallel$acc[(burnin_rounds+1):repeats,] / M)

  ## without random initialization
  ###############################################
  set.seed(i)
  ###
  priormean_parallel_without <- list()
  for (k in 1:nb.k) {
    priormean_parallel_without[[k]] <- rep(q/p, p)
  }
  L_parallel_without <- list()
  for (k in 1:nb.k) {
    L_parallel_without[[k]] <- p
  }
  
  start.time <- Sys.time()
  output_parallel_without = MAdaSub_parallel(data=data,nb.k=nb.k,M=M,repeats=repeats,priormean=priormean_parallel_without,L=L_parallel_without,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=FALSE,a_prior=NULL,b_prior=NULL,
                                     frac.par.update=1, burnin_rounds=burnin_rounds, numCores = numCores) 
  end.time <- Sys.time()
  Times_MAdaSub_parallel_without[i] <- difftime(end.time, start.time, unit = "sec")
  Estimates_MAdaSub_parallel_without[i,] <-   apply(output_parallel_without$importance.individual[repeats+1,,], 2, mean)
  Accrate_MAdaSub_parallel_without[i] <- mean(output_parallel_without$acc[(burnin_rounds+1):repeats,] / M)
  
  
  ## only r_j random initialization
  ###############################################
  set.seed(i)
  ###
  priormean_parallel_without <- list()
  for (k in 1:nb.k) {
    priormean_parallel[[k]] <- rep(runif(1,2,10)/p, p)
  }
  L_parallel_without <- list()
  for (k in 1:nb.k) {
    L_parallel_without[[k]] <- p
  }
  
  start.time <- Sys.time()
  output_parallel_only_rj = MAdaSub_parallel(data=data,nb.k=nb.k,M=M,repeats=repeats,priormean=priormean_parallel,L=L_parallel_without,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=FALSE,a_prior=NULL,b_prior=NULL,
                                             frac.par.update=1, burnin_rounds=burnin_rounds, numCores = numCores) 
  end.time <- Sys.time()
  Times_MAdaSub_parallel_only_rj[i] <- difftime(end.time, start.time, unit = "sec")
  Estimates_MAdaSub_parallel_only_rj[i,] <-   apply(output_parallel_only_rj$importance.individual[repeats+1,,], 2, mean)
  Accrate_MAdaSub_parallel_only_rj[i] <- mean(output_parallel_only_rj$acc[(burnin_rounds+1):repeats,] / M)
  
  
  ####################################################
  set.seed(i)
  start.time <- Sys.time()
  #outputMC3=MC3(data=data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,
  #              plot=FALSE, swap = swap)
  outputMC3 = mclapply(rep(Iterations, nb.k), MC3, mc.cores = numCores, data=data, priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,
                         plot=FALSE, swap = swap, S.start = NULL, burnin = burnin, 
                         mc.set.seed = TRUE,  mc.preschedule = TRUE )
  end.time <- Sys.time()
  MC3_estimates <- matrix(NA, nb.k, p)
  MC3_accrates <- rep(NA, nb.k)
  for (k in 1:nb.k) {
    MC3_estimates[k,] <- outputMC3[[k]]$importance.final
    MC3_accrates[k] <- outputMC3[[k]]$acc.prob.after.burnin
  }
  Times_MC3[i] <- difftime(end.time, start.time, unit = "sec")
  Estimates_MC3[i,] <- colMeans(MC3_estimates)
  Accrate_MC3[i] <- mean(MC3_accrates)
  
  #print(Estimates_MAdaSub[i,1:20])
  print("First 20 estimates MAdaSub parallel")
  print(Estimates_MAdaSub_parallel[i,1:20])
  print("First 20 estimates MAdaSub parallel without random init")
  print(Estimates_MAdaSub_parallel_without[i,1:20])
  print("First 20 estimates MAdaSub parallel only r_j random")
  print(Estimates_MAdaSub_parallel_only_rj[i,1:20])
  print("First 20 estimates MC3")
  print(Estimates_MC3[i,1:20])
  print(paste("Acceptance rate MAdaSub parallel:", Accrate_MAdaSub_parallel[i]))
  print(paste("Acceptance rate MAdaSub parallel without:", Accrate_MAdaSub_parallel_without[i]))
  print(paste("Acceptance rate MAdaSub parallel only r_j random:", Accrate_MAdaSub_parallel_only_rj[i]))
  print(paste("Acceptance rate MC3:", Accrate_MC3[i]))
  print("#######################################")
}

results = list(Estimates_MAdaSub_parallel_without = Estimates_MAdaSub_parallel_without,
               Estimates_MAdaSub_parallel = Estimates_MAdaSub_parallel,
               Estimates_MAdaSub_parallel_only_rj = Estimates_MAdaSub_parallel_only_rj,
               Estimates_MC3 = Estimates_MC3,
               Accrate_MAdaSub_parallel_without = Accrate_MAdaSub_parallel_without,
               Accrate_MAdaSub_parallel = Accrate_MAdaSub_parallel,
               Accrate_MAdaSub_parallel_only_rj = Accrate_MAdaSub_parallel_only_rj,
               Accrate_MC3 = Accrate_MC3,
               Times_MAdaSub_parallel_without = Times_MAdaSub_parallel_without,
               Times_MAdaSub_parallel = Times_MAdaSub_parallel,
               Times_MAdaSub_parallel_only_rj = Times_MAdaSub_parallel_only_rj,
               Times_MC3 = Times_MC3,
               n = n,
               p = p,
               SNR = SNR,
               M = M,
               burnin_rounds = burnin_rounds,
               repeats = repeats,
               Iterations = Iterations,
               burnin = burnin,
               numCores = numCores,
               nb.k = nb.k)

save(results, file=paste("MAdaSub_results_highdim_parallel_without2_sim_n", n, "_p", p, "_SNR", SNR, "_CPU", nb.k, "_Rep", repeats, "_M", M, ".RData", sep=""))



