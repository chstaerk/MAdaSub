
library("BAS")
library("coda")

RNGversion("4.0.5")

#setwd("C:\\Users\\Staerk\\scieboBonn\\MAdaSub PAPER\\R files")
source("MAdaSub_LOAD.R")

#devtools::install_github("albcab/scaleBVS")
library("scaleBVS")

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
  } else {
    sw <- rep(0, Iter)
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
  corr=0.9 
}

if (mychoice==3){
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

n_sim = 200
Iterations=20000
precision = 0.05

Ljs <- c(p/20, p/10, p/5, p/2, p, 2*p, 5*p, 10*p, 20*p)

acceptance_rates_MAdaSub = matrix(NA, n_sim, length(Ljs))
Iterations_converged_MAdaSub = matrix(NA, n_sim, length(Ljs))

#set.seed(22)
start.time <- Sys.time()
for (k in 1:n_sim) {
  #s0 = 5
  set.seed(k*22)
  print(paste("Simulated dataset", k))
  s0 = sample(0:10,1)
  beta1=numeric(p)
  S0=sort (sample(1:p,s0))
  beta1[S0]= runif(s0,-2,2)
  
  data=simdata(n,p,beta=beta1,corr= corr,family=family,sigma.normal=sigma.normal,blocks=blocks, corr.type=corr.type) 
  
  dataf = as.data.frame(data)
  out = bas.lm(y ~ .,
               data=dataf, prior="g-prior", #g-prior 
               modelprior=Bernoulli(priorprob),
               method="BAS", n.models = 2^p)
  
  # MAdaSub
  for (l in 1:length(Ljs)) {
    output=MAdaSub(data,Iter=Iterations,priormean=priormean, L=Ljs[l],const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE,automatic.stop=FALSE)
    
    acceptance_rates_MAdaSub[k,l] = output$acc/Iterations
    
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
    Iterations_converged_MAdaSub[k,l] = Iter_converged_MAdaSub
  }
  
  print(paste("Acceptance rate MAdaSub L=", Ljs , ":", acceptance_rates_MAdaSub[k,]))
  print(paste("Iterations for convergence MAdaSub L=", Ljs , ":", Iterations_converged_MAdaSub[k,]))
  print("#######################################")
}

end.time <- Sys.time()
end.time - start.time

if (corr==0.3) {
    results_lowdim_03 = list(acceptance_rates_MAdaSub = acceptance_rates_MAdaSub,
                             Iterations_converged_MAdaSub = Iterations_converged_MAdaSub)
    results = results_lowdim_03
  }

if (corr==0.9) {
    results_lowdim_09 = list(acceptance_rates_MAdaSub = acceptance_rates_MAdaSub,
                             Iterations_converged_MAdaSub = Iterations_converged_MAdaSub)
    results = results_lowdim_09
}

if (corr==0.99) {
  results_lowdim_099 = list(acceptance_rates_MAdaSub = acceptance_rates_MAdaSub,
                            Iterations_converged_MAdaSub = Iterations_converged_MAdaSub)
  results = results_lowdim_099
}


Iterations = 20000
p = 20
Ljs = c(p/20, p/10, p/5, p/2, p, 2*p, 5*p, 10*p, 20*p)

par(mfcol=c(2,3))
boxplot(results_lowdim_03_Lj$acceptance_rates_MAdaSub, names=Ljs, ylim=c(0,1))
boxplot(results_lowdim_03_Lj$Iterations_converged_MAdaSub, names=Ljs, ylim = c(0,Iterations+1))

boxplot(results_lowdim_09_Lj$acceptance_rates_MAdaSub, names=Ljs, ylim=c(0,1))
boxplot(results_lowdim_09_Lj$Iterations_converged_MAdaSub, names=Ljs, ylim = c(0,Iterations+1))

boxplot(results_lowdim_099_Lj$acceptance_rates_MAdaSub, names=Ljs, ylim=c(0,1))
boxplot(results_lowdim_099_Lj$Iterations_converged_MAdaSub, names=Ljs, ylim = c(0,Iterations+1))

# SAVE
#save(results_lowdim_03_Lj, file="results_lowdim_03_Lj.RData")
#save(results_lowdim_09_Lj, file="results_lowdim_09_Lj.RData")
#save(results_lowdim_099_Lj, file="results_lowdim_099_Lj.RData")

#load("results_lowdim_03_Lj.RData")
#load("results_lowdim_09_Lj.RData")
#load("results_lowdim_099_Lj.RData")

methodnames = c("MAdaSub (a)", "MAdaSub (b)", "MC3", "MC3 swap", "wTGS")


font = 1
par(las=1,lwd=1,cex.main=font, cex.lab=font, cex.axis=1.3)
par(cex.main=1.2)
options(scipen=5)

dev.off()
win.graph(width = 14, height = 8, pointsize = 8)

par(mfcol=c(2,3))
#par(mfcol=c(3,1))
par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))
results = results_lowdim_03_Lj
boxplot(results$acceptance_rates_MAdaSub, ylim = c(0,1), names = methodnames)
#expression(paste("Acceptance rates", rho," = 0.3"))
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, results$acceptance_rates_BVS, main = expression(bold(atop("Acceptance rate", rho ~" = 0.3"))), ylim = c(0,1), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3,  results$Iterations_converged_BVS, results$Iterations_converged_RB_BVS, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.3"))), ylim = c(0,Iterations+1), names = c(methodnames, "RB estimate"))
#boxplot(results$effective_sizes_median_MAdaSub, results$effective_sizes_median_PO, results$effective_sizes_median_MC3,  results$effective_sizes_median_BVS, main = expression(bold(atop("Median effective sample size", rho ~" = 0.3"))), ylim = c(0,Iterations+1), names = methodnames)

results = results_lowdim_09_Lj
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, results$acceptance_rates_MC3_swap, results$acceptance_rates_BVS, main = expression(bold(atop("Acceptance rate", rho ~" = 0.9"))), ylim = c(0,1), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, results$Iterations_converged_MC3_swap, results$Iterations_converged_BVS, results$Iterations_converged_RB_BVS, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.9"))), ylim = c(0,Iterations+1), names = c(methodnames, "RB estimate"))

results = results_lowdim_099_Lj
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, results$acceptance_rates_BVS, main = expression(bold(atop("Acceptance rate", rho ~" = 0.99"))), ylim = c(0,1), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3,  results$Iterations_converged_BVS, results$Iterations_converged_RB_BVS, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.99"))), ylim = c(0,Iterations+1), names = c(methodnames, "RB estimate"))

par(mar=c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

