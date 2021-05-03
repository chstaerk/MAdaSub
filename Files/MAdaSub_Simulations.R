
library("BAS")
library("coda")

#setwd("C:\\Users\\Staerk\\scieboBonn\\MAdaSub PAPER\\R files")
#source("MAdaSub_LOAD.R")

MC3 <-function (data,Iter,priormean,L,const=0,savings=1,family="normal",epsilon=1e-06,priorprob,prior,g,hyper=FALSE,a_prior=NULL,b_prior=NULL,plot=TRUE) { 
  
  p=ncol(data$x)
  n=nrow(data$x)
  
  CV=matrix(numeric(Iter/savings*p),p,floor(Iter/savings))
  CS=matrix(numeric(Iter/savings*p),p,floor(Iter/savings))
  
  CV.cur = numeric(p)
  CS.cur = numeric(p)
  
  V.size<-numeric(Iter)
  S.size<-numeric(Iter)
  
  values <- rep(NA,Iter)
  
  best.models = vector("list", Iter)
  
  S = integer(0)
  
  acc.prob = numeric(Iter)
  acc = 0
  
  for (t in 1:Iter) {
    
    j = sample(1:p,1)
    V = S
    if (is.element(j,S)) {
      V = setdiff(V, j)
    } else {
      V = union(V, j)
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
    
    if (values[t]==min(values[1:t])) best.S = S
    
    best.models[[t]] = S
    
    acc.prob[t] = acc/t
    
    if ((t %% 10000 == 0) && plot){ 
      #print(t)
      par(mfrow=c(1,1))
      plot(values,pch=20,main="MAdaSub" ,xlab="Iter",ylab="value")
      abline(h=EBIC_glm(data,which(relfreq>0.5),const,family=family),col="blue")
      abline(h=EBIC_glm(data,which(relfreq>0.9),const,family=family),col="red")
    }
  }
  
  best.S = best.models[[which.min(values)]]
  
  helpmat = matrix(rep((1:floor(Iter/savings))*savings,p) , p, floor(Iter/savings),byrow=TRUE)
  
  importance.hist = CS/helpmat 
  importance.final = CS.cur/Iter
  
  par(mfrow=c(1,1))
  if (plot){
    plot(values,pch=20,main="MAdaSub" ,xlab="Iter",ylab="value")
  }
  
  return(list(importance.final=importance.final,importance.hist=importance.hist,S.size=S.size,V.size=V.size,values=values,best.S=best.S,best.models=best.models, acc.prob= acc.prob, acc=acc, CS=CS))
}


family="normal"
n=60 # 60
p=20 #20

sigma.normal= 1 #1
corr.type="toeplitz"

# Specify correlation between adjacent covariates in Toeplitz correlation 
mychoice <- menu( c("rho=0.3                                  ", 
                    "rho=0.9                                  ", 
                    "rho=0.99                                 "), 
                    graphics=TRUE, title="Choose Toeplitz correlation " )

if (mychoice==1){ 
  corr=0.3 
} 

if (mychoice==2){
  corr.type="global"
  corr=0.9 
}

if (mychoice==3){
  corr.type="toeplitz"
  corr=0.99 
}

#corr=0.9  #0.3 0.9 0.99

##################################################################.
# Begin MAdaSub algo
##################################################################

const=1
savings = 1 #1
#K=1 # n

q=10

epsilon = 1/p #1e-6
priormean = q/p #1/2 # q/p
L = p #rep(p,p)  #K=p/L
priorprob= q/p #1/2
g=n

prior="gprior" # g-prior (original form with Jeffreys prior on variance)
#prior="EBIC"
#prior="independence" 
hyper = FALSE
a_prior = NULL
b_prior = NULL

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)


choices <-c( "100 simulated datasets (as in paper)", "Custom number of simulated datasets")
mychoice <- menu( choices , graphics=TRUE, title="How many simulated datasets should be considered?" )

if (mychoice==1) {
  n_sim=100 
}
if (mychoice==2) {
  n_sim = -1
  while(n_sim < 1 ){
    n_sim <- readline("Enter the number of simulated datasets: ")
    n_sim <- ifelse(grepl("\\D",n_sim),-1,as.integer(n_sim))
    if(is.na(n_sim)){break}  # breaks when hit enter
  }
}

#n_sim = 100
Iterations=20000
precision = 0.05

effective_sizes_MAdaSub = matrix(NA,nrow = n_sim, ncol = p)
effective_sizes_PO = matrix(NA,nrow = n_sim, ncol = p)
effective_sizes_MC3 = matrix(NA,nrow = n_sim, ncol = p)

effective_sizes_median_MAdaSub = rep(NA, n_sim)
effective_sizes_median_PO = rep(NA, n_sim)
effective_sizes_median_MC3 = rep(NA, n_sim)

acceptance_rates_MAdaSub = rep(NA, n_sim)
acceptance_rates_PO = rep(NA, n_sim)
acceptance_rates_MC3 = rep(NA, n_sim)

Iterations_converged_MAdaSub = rep(NA, n_sim)
Iterations_converged_PO = rep(NA, n_sim)
Iterations_converged_MC3 = rep(NA, n_sim)

set.seed(22)
start.time <- Sys.time()
for (k in 1:n_sim) {
  #s0 = 5
  print(paste("Simulated dataset", k))
  s0 = sample(0:10,1)
  beta1=numeric(p)
  S0=sort (sample(1:p,s0))
  beta1[S0]= runif(s0,-2,2)
  
  data=simdata(n,p,beta=beta1,corr= corr,family=family,sigma.normal=sigma.normal,blocks=blocks, corr.type=corr.type) 
  
  output=MAdaSub(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE,automatic.stop=FALSE)
  
  outputMC3=MC3(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE)
  
  priormean2 = rep(q/p,p)
  PO = numeric(p)
  if (prior == "EBIC") {
    value_null_model = EBIC_glm(data, indices = c(), const=const,family=family ) 
  } else {
    value_null_model = marginal_log_like(data=data,indices=c(),family=family,prior=prior,g=g) + 
      log_modelprior(indices=c(),p=p,priorprob=priorprob,hyper=hyper,a_prior=a_prior,b_prior=b_prior) 
  }
  for (j in 1:p){
    if (prior == "EBIC") {
      PO[j] = -1/2*EBIC_glm(data, indices = j, const=const,family=family ) + 1/2*value_null_model
    } else {
      PO[j] = marginal_log_like(data=data,indices=j,family=family,prior=prior,g=g) + 
        log_modelprior(indices=j,p=p,priorprob=priorprob,hyper=hyper,a_prior=a_prior,b_prior=b_prior)  - 
        value_null_model
    }
    priormean2[j] = exp(PO[j])/(1+exp(PO[j]))
  } 
  priormean2[priormean2>0.9] = 0.9
  priormean2[priormean2<1/p] = 1/p
  
  output_PO=MAdaSub(data,Iter=Iterations,priormean=priormean2, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE,automatic.stop=FALSE)
  
  ##### Compute effective sample size #############################
  
  Gamma_indi_mat = output$CS
  for (i in (2:(dim(output$CS)[2]))) {
    Gamma_indi_mat[,i] = output$CS[,i] - output$CS[,i-1]
  }
  Gamma_indi_mcnc= as.mcmc(t(Gamma_indi_mat)) 
  eff_size_MAdaSub = effectiveSize(Gamma_indi_mcnc)
  
  Gamma_indi_mat = outputMC3$CS
  for (i in (2:(dim(outputMC3$CS)[2]))) {
    Gamma_indi_mat[,i] = outputMC3$CS[,i] - outputMC3$CS[,i-1]
  }
  Gamma_indi_mc3= as.mcmc(t(Gamma_indi_mat)) 
  eff_size_MC3 = effectiveSize(Gamma_indi_mc3)
  
  Gamma_indi_mat = output_PO$CS
  for (i in (2:(dim(output_PO$CS)[2]))) {
    Gamma_indi_mat[,i] = output_PO$CS[,i] - output_PO$CS[,i-1]
  }
  Gamma_indi_mc= as.mcmc(t(Gamma_indi_mat)) 
  eff_size_PO = effectiveSize(Gamma_indi_mc)

  effective_sizes_MAdaSub[k,] = eff_size_MAdaSub
  effective_sizes_median_MAdaSub[k] = median(eff_size_MAdaSub)
  
  effective_sizes_MC3[k,] = eff_size_MC3
  effective_sizes_median_MC3[k] = median(eff_size_MC3)
  
  effective_sizes_PO[k,] = eff_size_PO
  effective_sizes_median_PO[k] = median(eff_size_PO)
  
  acceptance_rates_MAdaSub[k] = output$acc/Iterations
  acceptance_rates_PO[k] = output_PO$acc/Iterations
  acceptance_rates_MC3[k] = outputMC3$acc/Iterations
  
  ###########################
    dataf = as.data.frame(data)
    out = bas.lm(y ~ .,
                 data=dataf, prior="g-prior", #g-prior 
                 modelprior=Bernoulli(priorprob),
                 method="BAS", n.models = 2^p)
    
    #outRS = bas.lm(y ~ .,
    #             data=dataf, prior="g-prior", #g-prior 
    #             modelprior=Bernoulli(priorprob),
    #             method="MCMC", n.models = Iterations)
    
  
  Iter_adj = length(output$importance.hist[1,])
  converged = logical(Iter_adj)
  for (t in 1:Iter_adj){
    if (all(abs(output$importance.hist[,t] - out$probne0[-1])<= precision)) converged[t] = TRUE
  }
  if (all(!converged)) {
    Iter_converged_MAdaSub = Iterations + 1 
  } else {
    Iter_converged_MAdaSub = which.max(converged)
  }
    
  Iter_adj = length(output_PO$importance.hist[1,])
  converged = logical(Iter_adj)
  for (t in 1:Iter_adj){
    if (all(abs(output_PO$importance.hist[,t] - out$probne0[-1])<= precision)) converged[t] = TRUE
  }
  if (all(!converged)) {
    Iter_converged_PO = Iterations + 1 
  } else {
    Iter_converged_PO = which.max(converged)
  }
  
  Iter_adj = length(outputMC3$importance.hist[1,])
  converged = logical(Iter_adj)
  for (t in 1:Iter_adj){
    if (all(abs(outputMC3$importance.hist[,t] - out$probne0[-1])<= precision)) converged[t] = TRUE
  }
  if (all(!converged)) {
    Iter_converged_MC3 = Iterations + 1 
  } else {
    Iter_converged_MC3 = which.max(converged)
  }
  
  Iterations_converged_MAdaSub[k] = Iter_converged_MAdaSub
  Iterations_converged_PO[k] = Iter_converged_PO
  Iterations_converged_MC3[k] = Iter_converged_MC3
  print(paste("Acceptance rate MAdaSub (a):", output$acc/Iterations))
  print(paste("Acceptance rate MAdaSub (b):", output_PO$acc/Iterations))
  print(paste("Acceptance rate MC3:", outputMC3$acc/Iterations))
  print(paste("Iterations for convergence MAdaSub (a):", Iter_converged_MAdaSub))
  print(paste("Iterations for convergence MAdaSub (b):", Iter_converged_PO))
  print(paste("Iterations for convergence MC3:", Iter_converged_MC3))
  print("#######################################")
}

end.time <- Sys.time()
end.time - start.time

if (corr==0.3) {
    results_lowdim_03 = list(acceptance_rates_MAdaSub = acceptance_rates_MAdaSub,
                        acceptance_rates_MC3 = acceptance_rates_MC3,
                        acceptance_rates_PO = acceptance_rates_PO,
                        effective_sizes_median_MAdaSub = effective_sizes_median_MAdaSub,
                        effective_sizes_median_PO = effective_sizes_median_PO,
                        effective_sizes_median_MC3 = effective_sizes_median_MC3,
                        Iterations_converged_MAdaSub = Iterations_converged_MAdaSub,
                        Iterations_converged_PO = Iterations_converged_PO,
                        Iterations_converged_MC3 = Iterations_converged_MC3)
    results = results_lowdim_03
  }

if (corr==0.9) {
    results_lowdim_09 = list(acceptance_rates_MAdaSub = acceptance_rates_MAdaSub,
                          acceptance_rates_MC3 = acceptance_rates_MC3,
                          acceptance_rates_PO = acceptance_rates_PO,
                          effective_sizes_median_MAdaSub = effective_sizes_median_MAdaSub,
                          effective_sizes_median_PO = effective_sizes_median_PO,
                          effective_sizes_median_MC3 = effective_sizes_median_MC3,
                          Iterations_converged_MAdaSub = Iterations_converged_MAdaSub,
                          Iterations_converged_PO = Iterations_converged_PO,
                          Iterations_converged_MC3 = Iterations_converged_MC3  )
    results = results_lowdim_09
}

if (corr==0.99) {
  results_lowdim_099 = list(acceptance_rates_MAdaSub = acceptance_rates_MAdaSub,
                           acceptance_rates_MC3 = acceptance_rates_MC3,
                           acceptance_rates_PO = acceptance_rates_PO,
                           effective_sizes_median_MAdaSub = effective_sizes_median_MAdaSub,
                           effective_sizes_median_PO = effective_sizes_median_PO,
                           effective_sizes_median_MC3 = effective_sizes_median_MC3,
                           Iterations_converged_MAdaSub = Iterations_converged_MAdaSub,
                           Iterations_converged_PO = Iterations_converged_PO,
                           Iterations_converged_MC3 = Iterations_converged_MC3  )
  results = results_lowdim_099
}


# SAVE
#save(results_lowdim_03, file="results_lowdim_03_epsilon.RData")
#save(results_lowdim_09, file="results_lowdim_09_epsilon.RData")
#save(results_lowdim_099, file="results_lowdim_099_epsilon.RData")

#load("results_lowdim_03_epsilon.RData")
#load("results_lowdim_09_epsilon.RData")
#load("results_lowdim_099_epsilon.RData")

methodnames = c("MAdaSub (a)", "MAdaSub (b)", "MC3")


font = 1
par(las=1,lwd=1,cex.main=font, cex.lab=font, cex.axis=1.3)
par(cex.main=1.2)
options(scipen=5)

dev.off()
win.graph(width = 14, height = 8, pointsize = 8)

#par(mfcol=c(2,3))
par(mfcol=c(2,1))
par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))
#results = results_lowdim_03
#expression(paste("Acceptance rates", rho," = 0.3"))
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, main = expression(bold(atop("Acceptance rate", rho ~" = 0.3"))), ylim = c(0,1), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.3"))), ylim = c(0,Iterations+1), names = methodnames)

#results = results_lowdim_09
#boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, main = expression(bold(atop("Acceptance rate", rho ~" = 0.9"))), ylim = c(0,1), names = methodnames)
#boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.9"))), ylim = c(0,Iterations+1), names = methodnames)

#results = results_lowdim_099
#boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, main = expression(bold(atop("Acceptance rate", rho ~" = 0.99"))), ylim = c(0,1), names = methodnames)
#boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.99"))), ylim = c(0,Iterations+1), names = methodnames)
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

