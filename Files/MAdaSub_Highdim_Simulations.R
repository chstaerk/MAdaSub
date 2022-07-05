
library("BAS")
library("coda")

library("glmnet")

#RNGkind(kind = "default", normal.kind = "default", sample.kind = "default")
RNGkind("L'Ecuyer-CMRG")

#setwd("C:\\Users\\Staerk\\scieboBonn\\MAdaSub PAPER\\R files")
#source("MAdaSub_LOAD.R")

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

# Specify sample size n 
mychoice <- menu( c("n = 500                                 ", 
                    "n = 1000                                 "), 
                  graphics=TRUE, title="Choose sample size " )
if (mychoice==1){ 
  n = 500  
} 
if (mychoice==2){
  n = 1000 
}

# Specify number of covariates
mychoice <- menu( c("p = 500                                 ", 
                    "p = 5000                                 "), 
                  graphics=TRUE, title="Choose number of covariates " )
if (mychoice==1){ 
  p = 500  
} 
if (mychoice==2){
  p = 5000 
}


# Specify SNR
mychoice <- menu( c("SNR = 0.5                                 ", 
                    "SNR = 1                                 ",
                    "SNR = 2                                 ", 
                    "SNR = 3                                 "), 
                  graphics=TRUE, title="Choose signal-to-noise ratio " )
if (mychoice==1){ 
  SNR = 0.5 
} 
if (mychoice==2){
  SNR = 1 
}
if (mychoice==3){ 
  SNR = 2 
} 
if (mychoice==4){
  SNR = 3 
}


choices <-c( "200 runs (as in paper)", "Custom number of runs")
mychoice <- menu( choices , graphics=TRUE, title="How many runs of the algorithms should be considered?" )

if (mychoice==1) {
  N_sim=200 
}
if (mychoice==2) {
  N_sim = -1
  while(N_sim < 1 ){
    N_sim <- readline("Enter the number of runs: ")
    N_sim <- ifelse(grepl("\\D",N_sim),-1,as.integer(N_sim))
    if(is.na(N_sim)){break}  # breaks when hit enter
  }
}

choices <-c( "5 CPUs (as in paper)", "1 CPU (serial computations)" , "Custom number of CPUs")
mychoice <- menu( choices , graphics=TRUE, title="How many CPUs should be used?" )

if (mychoice==1) {
  numCores <- 5  
}
if (mychoice==2) {
  numCores <- 1 
}
if (mychoice==3) {
  numCores = -1
  while(numCores < 1 ){
    numCores <- readline("Enter the number of CPUs: ")
    numCores <- ifelse(grepl("\\D",numCores),-1,as.integer(numCores))
    if(is.na(numCores)){break}  # breaks when hit enter
  }
}

choices <-c( "5 chains (as in paper)", "Custom number of chains (>=2)")
mychoice <- menu( choices , graphics=TRUE, title="How many parallel chains for each algorithm should be used?" )

if (mychoice==1) {
  nb.k <- 5  
}
if (mychoice==2) {
  nb.k = -1
  while(nb.k < 1 ){
    nb.k <- readline("Enter the number of chains: ")
    nb.k <- ifelse(grepl("\\D",nb.k),-1,as.integer(nb.k))
    if(is.na(nb.k)){break}  # breaks when hit enter
  }
}


#n = 500 
#p = 5000
#SNR = 0.5  

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


##################################################################
# Begin MAdaSub algo
##################################################################

q=10 

epsilon = 1/p
priormean = q/p 
L = p 
priorprob= q/p 
g=9 

if (p==500) M <- 1000 
if (p==5000) M <- 10000 

burnin_rounds <- 25 # 25
repeats <- 50 # 50 
Iterations <- M * repeats  
burnin <- M * burnin_rounds 
#nb.k <- 5 
savings = 100
#numCores <- 5 # nb.k

prior="independence" 
hyper = FALSE
a_prior = NULL
b_prior = NULL

swap = TRUE

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)

#N_sim <- 200 

Times_MAdaSub <- rep(NA,N_sim)
Times_MAdaSub_parallel <- rep(NA,N_sim)
Times_MAdaSub_serial <- rep(NA,N_sim)
Times_MC3 <- rep(NA,N_sim)

Estimates_MAdaSub <- matrix(NA, N_sim, p)
Estimates_MAdaSub_parallel <- matrix(NA, N_sim, p)
Estimates_MAdaSub_serial <- matrix(NA, N_sim, p)
Estimates_MC3 <- matrix(NA, N_sim, p)

Accrate_MAdaSub <- rep(NA,N_sim)
Accrate_MAdaSub_parallel <- rep(NA,N_sim)
Accrate_MAdaSub_serial <- rep(NA,N_sim)
Accrate_MC3 <- rep(NA,N_sim)

for (i in 1:N_sim){
  print(paste("Run ", i))
  set.seed(i)
  
  priormean_serial <- q/p
  start.time <- Sys.time()
  #output=MAdaSub(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,
  #               plot=FALSE,automatic.stop=FALSE, burnin=burnin)
  output = mclapply(rep(Iterations, nb.k), MAdaSub, mc.cores = numCores, data=data, priormean=priormean_serial, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,
                    plot=FALSE,automatic.stop=FALSE, burnin=burnin,
                    mc.set.seed = TRUE,  mc.preschedule = TRUE )
  end.time <- Sys.time()
  MAdaSub_estimates <- matrix(NA, nb.k, p)
  MAdaSub_accrates <- rep(NA, nb.k)
  for (k in 1:nb.k) {
    MAdaSub_estimates[k,] <- output[[k]]$importance.final
    MAdaSub_accrates[k] <- output[[k]]$acc.prob.after.burnin
  }
  Times_MAdaSub[i] <- difftime(end.time, start.time, unit = "sec")
  Estimates_MAdaSub[i,] <- colMeans(MAdaSub_estimates)
  Accrate_MAdaSub[i] <- mean(MAdaSub_accrates)
  
  
  ###############################################
  set.seed(i)
  ###
  priormean_parallel <- list()
  for (k in 1:nb.k) {
    #priormean_parallel[[k]] <- rep(q/p, p)
    priormean_parallel[[k]] <- rep(runif(1,2,10)/p, p)
  }
  L_parallel <- list()
  for (k in 1:nb.k) {
    #L_parallel[[k]] <- p
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
  
  ####################################################
  set.seed(i)
  start.time <- Sys.time()
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
  print("First 20 PIP estimates MAdaSub parallel")
  print(round(Estimates_MAdaSub_parallel[i,1:20], 6))
  print("First 20 PIP estimates MAdaSub serial")
  print(round(Estimates_MAdaSub[i,1:20], 6))
  print("First 20 PIP estimates MC3")
  print(round(Estimates_MC3[i,1:20], 6))
  print(paste("Acceptance rate MAdaSub parallel: ", round(100*Accrate_MAdaSub_parallel[i], 2), "%", sep=""))
  print(paste("Acceptance rate MAdaSub serial: ", round(100*Accrate_MAdaSub[i], 2), "%", sep=""))
  print(paste("Acceptance rate MC3: ", round(100*Accrate_MC3[i], 2), "%", sep=""))
  print("#######################################")
}

results = list(Estimates_MAdaSub = Estimates_MAdaSub,
               Estimates_MAdaSub_parallel = Estimates_MAdaSub_parallel,
               Estimates_MC3 = Estimates_MC3,
               Accrate_MAdaSub = Accrate_MAdaSub,
               Accrate_MAdaSub_parallel = Accrate_MAdaSub_parallel,
               Accrate_MC3 = Accrate_MC3,
               Times_MAdaSub = Times_MAdaSub,
               Times_MAdaSub_parallel = Times_MAdaSub_parallel,
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

save(results, file=paste("MAdaSub_results_highdim_sim_n", n, "_p", p, "_SNR", SNR, "_CPU", numCores, "_Rep", repeats, "_M", M, ".RData", sep=""))



var_MAdaSub <- apply(results$Estimates_MAdaSub, 2, var)
var_MAdaSub_parallel <- apply(results$Estimates_MAdaSub_parallel, 2, var)
var_MC3 <- apply(results$Estimates_MC3, 2, var)

r_hat <- var_MC3 * median(results$Times_MC3) / (var_MAdaSub * median(results$Times_MAdaSub))
r_hat_parallel <- var_MC3 * median(results$Times_MC3) / (var_MAdaSub_parallel * median(results$Times_MAdaSub_parallel))

n_best <- 20
highest_PIP <- order(colMeans(rbind(results$Estimates_MAdaSub_parallel,
                                    results$Estimates_MAdaSub,
                                    results$Estimates_MC3)), decreasing=TRUE)[1:n_best]

sorted_PIP <- sort(colMeans(rbind(results$Estimates_MAdaSub_parallel,
                                  results$Estimates_MAdaSub,
                                  results$Estimates_MC3)), decreasing=TRUE)

print("#######################################")
print("Results:")
print("#######################################")
print(paste("Median acceptance rate MAdaSub: ", round(median(results$Accrate_MAdaSub)*100, 2), "%", sep=""))
print(paste("Median acceptance rate MAdaSub parallel: ", round(median(results$Accrate_MAdaSub_parallel)*100, 2), "%", sep=""))
print(paste("Median acceptance rate MC3: ", round(median(results$Accrate_MC3)*100, 2), "%", sep=""))
print("#######################################")
print(paste("Median r hat (top 20 variables) MAdaSub vs. MC3: ", round(median(r_hat[highest_PIP], na.rm=TRUE),1)))
print(paste("Median r hat (top 20 variables) MAdaSub parallel vs. MC3: ", round(median(r_hat_parallel[highest_PIP], na.rm=TRUE),1)))
print("#######################################")
print(paste("Median r hat (all variables) MAdaSub vs. MC3: ", round(median(r_hat, na.rm=TRUE),1)))
print(paste("Median r hat (all variables) MAdaSub parallel vs. MC3: ", round(median(r_hat_parallel, na.rm=TRUE),1)))

