
library("BAS")
library("coda")

#RNGkind(kind = "default", normal.kind = "default", sample.kind = "default")
RNGkind("L'Ecuyer-CMRG")

#setwd("C:\\Users\\Staerk\\scieboBonn\\MAdaSub PAPER\\R files")
#source("MAdaSub_LOAD.R")

IndMetr <-function (data,Iter,priormean,L,const=0,savings=1,family="normal",epsilon=1e-06,priorprob,prior,g,hyper=FALSE,a_prior=NULL,b_prior=NULL,plot=TRUE) { 
  
  p=ncol(data$x)
  n=nrow(data$x)
  
  relfreq=numeric(p)
  relfreq[1:p] = priormean
  
  relfreq[relfreq>1-epsilon] = 1-epsilon
  relfreq[relfreq<epsilon] = epsilon
  
  CV=matrix(numeric(Iter/savings*p),p,floor(Iter/savings))
  CS=matrix(numeric(Iter/savings*p),p,floor(Iter/savings))
  
  CV.cur = numeric(p)
  CS.cur = numeric(p)
  
  V.size<-numeric(Iter)
  S.size<-numeric(Iter)
  
  values <- rep(NA,Iter)
  
  best.models = vector("list", Iter)
  
  b = rbinom(p,1,relfreq)
  S = which(b==1)
  #S = integer(0)
  
  acc.prob = numeric(Iter)
  acc = 0
  
  for (t in 1:Iter) {
    
    b = rbinom(p,1,relfreq)
    V = which(b==1)
    
    #if (length(V)>U_C){ V=sample(V,U_C) } 
    
    if (prior=="EBIC"){
      alpha = (-1/2)* EBIC_glm(data=data,indices=V,const=const,family=family) +
        log_prob(indices=S,p=p,relfreq=relfreq) - 
        (-1/2)* EBIC_glm(data=data,indices=S,const=const,family=family) -
        log_prob(indices=V,p=p,relfreq=relfreq) } else {
      alpha =  marginal_log_like(data=data,indices=V,family=family,prior=prior,g=g) + 
        log_modelprior(indices=V,p=p,priorprob=priorprob) +
        log_prob(indices=S,p=p,relfreq=relfreq) - 
        marginal_log_like(data=data,indices=S,family=family,prior=prior,g=g)  - 
        log_modelprior(indices=S,p=p,priorprob=priorprob) -
        log_prob(indices=V,p=p,relfreq=relfreq)
    }
    u=log(runif(1))
    
    if (u<=alpha) { 
      S = V 
      acc = acc + 1
    }
    
    if (prior=="EBIC"){
      values[t] = (-1/2) * EBIC_glm(data=data,indices=S,const=const,family=family) } else {
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
  importance = CS.cur/Iter
  
  helpmat = matrix(rep((1:floor(Iter/savings))*savings,p) , p, floor(Iter/savings),byrow=TRUE)
  
  relfreq.hist=(L*priormean + CS)/(L+helpmat) 
  relfreq.final = relfreq
  #relfreq.final=relfreq.hist[,floor(Iter/savings)]
  
  
  par(mfrow=c(1,1))
  if (plot){
  plot(values,pch=20,main="MAdaSub" ,xlab="Iter",ylab="value")
  #abline(h=EBIC_glm(data,which(relfreq>0.5),const,family=family),col="blue")
  #abline(h=EBIC_glm(data,which(relfreq>0.9),const,family=family),col="red")
  }
  
  #return(list(importance=importance,values=values, acc.prob= acc.prob))
  return(list(relfreq.hist=relfreq.hist, relfreq.final=relfreq.final, importance=importance,S.size=S.size,V.size=V.size,values=values,best.S=best.S,best.models=best.models, acc.prob= acc.prob, acc=acc, CS=CS))
}

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
      #abline(h=EBIC_glm(data,which(relfreq>0.5),const,family=family),col="blue")
      #abline(h=EBIC_glm(data,which(relfreq>0.9),const,family=family),col="red")
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
n=60 
p=20 

beta1=numeric(p)
S0=c(1,2,3,4,5)
beta1[c(1,2,3,4,5)]=c(0.4,0.8,1.2,1.6,2)
S0=5
corr=0.9  
sigma.normal= 1 
blocks=10
corr.type="toeplitz"
set.seed(22)
data=simdata(n,p,beta=beta1,corr= corr,family=family,sigma.normal=sigma.normal,blocks=blocks, corr.type=corr.type)


print("Please wait a minute until all algorithms have been completed...")

##################################################################
# Begin MAdaSub algo
##################################################################

savings = 1 

q=10

epsilon = 1/p 
priormean = q/p 
L = p 
priorprob= q/p 
g=n

prior="gprior" # g-prior (original form with Jeffreys prior on variance)
hyper = FALSE
a_prior = NULL
b_prior = NULL

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)

Iterations=20000 
set.seed(22)
start.time <- Sys.time()
output=MAdaSub(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE,automatic.stop=FALSE)
end.time <- Sys.time()
end.time - start.time

#===============================================
# Different choices of r^0 
#===============================================

priormean = rep(q/p,p)
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
  priormean[j] = exp(PO[j])/(1+exp(PO[j]))
} 
priormean[priormean>0.9] = 0.9
priormean[priormean<1/p] = 1/p


start.time <- Sys.time()
output_PO=MAdaSub(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE,automatic.stop=FALSE)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
output_ind_false_PO=IndMetr(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE)
end.time <- Sys.time()
end.time - start.time


#===============================================

priormean = rep(q/p,p)
start.time <- Sys.time()
outputMC3=MC3(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE)
end.time <- Sys.time()
end.time - start.time

# Run BAS with full enumeration 
if (p<=35) {
dataf = as.data.frame(data)
start.time <- Sys.time()
out = bas.lm(y ~ .,
             data=dataf, prior="g-prior", #g-prior 
             modelprior=Bernoulli(priorprob),
             #initprobs= "eplogp", 
             method="BAS",
             n.models = 2^p)
end.time <- Sys.time()
end.time - start.time
}

#===========================================================================================

dataf = as.data.frame(data)
start.time <- Sys.time()
outMCMC = bas.lm(y ~ .,
             data=dataf, prior="g-prior", #g-prior 
             modelprior=Bernoulli(priorprob),
             #initprobs= "eplogp", 
             method="MCMC", MCMC.iterations=Iterations)
end.time <- Sys.time()
end.time - start.time


#===========================================================================================



# Non-adaptive sampler with prior marginals
priormean = q/p #1/2 
start.time <- Sys.time()
output_ind_false=IndMetr(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE)
end.time <- Sys.time()
end.time - start.time


# Non-adaptive sampler with true posterior marginals (from BAS)
priormean = out$probne0[-1]
start.time <- Sys.time()
output_ind_true=IndMetr(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE)
end.time <- Sys.time()
end.time - start.time


# output_ind_true$acc/Iterations


# Two additional runs of MAdaSub with different choices of variance parameter L
priormean = q/p # 1/2 

L=p/n
#q=5
set.seed(22)
start.time <- Sys.time()
output2=MAdaSub(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE)
end.time <- Sys.time()
end.time - start.time

L=100*p
#q=2

start.time <- Sys.time()
output3=MAdaSub(data,Iter=Iterations,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,plot=FALSE)
end.time <- Sys.time()
end.time - start.time



#############################################
## Plots for paper
#############################################

dev.off()
win.graph(width = 10, height = 6, pointsize = 10)

font = 1
par(las=1,lwd=1,cex.main=font, cex.lab=1.8, cex.axis=1.5)
par(cex.main=2)
par(mar=c(5, 4, 4, 1), mgp = c(2.5, 1, 0))
par(las=1)
options(scipen=5)

plotIter = 1:5000
#### Evolution of sampled values  
#selected = which(output$relfreq.hist[,Iterations/savings]>0.8) #output$best.S
ra1 = range(c(output$values,output_ind_true$values),na.rm=TRUE)
ra2 = range(c(output$values,output_ind_true$values,output_ind_false$values),na.rm=TRUE)
par(mfrow=c(4,1))
plot(output_ind_false$values[plotIter],pch=20,xlab="Iteration",ylab="",ylim=ra1,main="Non-adaptive sampler with prior marginals")
plot(output$values[plotIter],pch=20,xlab="Iteration",ylab="",ylim=ra1,main="MAdaSub")
#plot(output_PO$values[plotIter],pch=20,xlab="Iteration",ylab="",ylim=ra1,main="MAdaSub")
plot(output_ind_true$values[plotIter],pch=20,xlab="Iteration",ylab="",ylim=ra1,main="Non-adaptive sampler with posterior marginals")
plot(outputMC3$values[plotIter],pch=20,xlab="Iteration",ylab="",ylim=ra1,main="MC3")
# plot(output_ind_false_PO$values,pch=20,xlab="Iteration",ylab="",ylim=ra2,main="Non-adaptive sampler with prior marginals")


### Plot of sampled sizes of V^(t) and f(V^(t))

win.graph(width = 10, height = 6, pointsize = 10)
par(mfrow=c(4,1))
par(las=1)
plot(output_ind_false$V.size[plotIter],pch=20,ylim=c(0,max(output_ind_false$V.size)),xlab="Iteration",ylab="Model size",main="Non-adaptive sampler with prior marginals")
points(output_ind_false$S.size[plotIter],col='red',pch=20,cex=0.7)

#plot(output_ind_false_PO$V.size,pch=20,ylim=c(0,max(output_ind_false$V.size)),xlab="Iteration",ylab="Model size",main="Non-adaptive sampler with prior marginals")
#points(output_ind_false_PO$S.size,col='red',pch=20,cex=0.7)

plot(output$V.size[plotIter],pch=20,ylim=c(0,max(output_ind_false$V.size)),xlab="Iteration",ylab="Model size",main="MAdaSub")
points(output$S.size[plotIter],col='red',pch=20,cex=0.7)

#plot(output_PO$V.size,pch=20,ylim=c(0,max(output_ind_false$V.size)),xlab="Iteration",ylab="Model size",main="MAdaSub")
#points(output_PO$S.size,col='red',pch=20,cex=0.7)

plot(output_ind_true$V.size[plotIter],pch=20,ylim=c(0,max(output_ind_false$V.size)),xlab="Iteration",ylab="Model size",main="Non-adaptive sampler with posterior marginals")
points(output_ind_true$S.size[plotIter],col='red',pch=20,cex=0.7)

plot(outputMC3$V.size[plotIter],pch=20,ylim=c(0,max(output_ind_false$V.size)),xlab="Iteration",ylab="Model size",main="MC3")
points(outputMC3$S.size[plotIter],col='red',pch=20,cex=0.7)

### Plot of acceptance rate
win.graph(width = 10, height = 6, pointsize = 10)
par(mfrow=c(1,1))
par(las=1,lwd=1,cex.main=font, cex.lab=1.6, cex.axis=1.4)
par(cex.main=1.6)
par(las=1)
par(mar=c(4, 4, 4, 1))
plot(output$acc.prob,pch=20,ylim=c(0,1),xlab="Iteration",ylab="", main="Acceptance rate",lwd=2, type="l")
#points(output_PO$acc.prob,pch=20,col="yellow",cex=0.5)
lines(output_ind_true$acc.prob,pch=20,col="red",lwd=2)
lines(output_ind_false$acc.prob,pch=20,col="blue",lwd=2)
lines(output2$acc.prob,pch=20,col="orange",lwd=2)
lines(output3$acc.prob,pch=20,col="purple",lwd=2)
lines(outputMC3$acc.prob,pch=20,col="gray",lwd=2)
abline(h=output_ind_true$acc/Iterations, col="red",lwd=1)
abline(h=output_ind_false$acc/Iterations, col="blue",lwd=1)
legend(1000, 1, 
       c("Non-adaptive (prior)", 
         "MAdaSub (L=p)",
         "Non-adaptive (posterior)", 
         "MAdaSub (L=p/n)",
         "MC3", 
         "MAdaSub (L=100p)"),
       col = c("blue", "black", "red", "orange", "gray", "purple"),
       lty=1,
       ncol=3,
       lwd=2,
       cex=1.3,
       pt.cex = 1)


#plot(output$relfreq.final, output_PO$relfreq.final)
#abline(a=0,b=1)


# Automatic check for convergence
precision = 0.005
Iter.adj = length(output$relfreq.hist[1,])
converged = logical(Iter.adj)
for (t in 1:Iter.adj){
  if (all(abs(output$relfreq.hist[,t] - output$importance.hist[,t] )<= precision)) converged[t] = TRUE
}
Iter.converged = which.max(converged)


#par(mfrow=c(1,1))
#plot(converged,type ="l")

plotIter2 = 1:5000

win.graph(width = 10, height = 6, pointsize = 10)

############################################################################
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
layout(mat = m, heights = c(0.3,0.3,0.3,0.1))
############################################################################

par(mar=c(4, 3, 3, 1), mgp = c(2.5, 1, 0))
par(las=1)


#### Evolution of r_j^(t) and estimated PIP (only every "savings" iteration)
#par(mfrow=c(3,3))
zehner=0
for (i in 1:9){
  plot(output$relfreq.hist[10*zehner+i,plotIter2],pch=20,xlab="Iteration",ylab="",main=paste("Variable ",10*zehner+i, sep=""),ylim=c(0,1),cex=0.3, type="l", lwd=2, lty=2)
  lines(output$importance.hist[10*zehner+i,plotIter2],pch=20,col="gray",cex=0.3, lwd=2, lty=1)
  lines(output$relfreq.hist[10*zehner+i,plotIter2], type="l", lwd=2, lty=2)
  abline(v = Iter.converged)
  #abline(h=out$probne0[-1][i],col="red") 
}

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend = c("Proposal probabilities", "Empirical inclusion frequencies"), 
       col=c("black", "gray"), cex=2, box.lty=1, 
       pt.bg = 'white', lty=c(2,1), lwd =2, ncol=2) #,
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))


#### Evolution of r_j^(t) (only every "savings" iteration)

win.graph(width = 10, height = 6, pointsize = 10)

############################################################################
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
layout(mat = m, heights = c(0.3,0.3,0.3,0.1))
############################################################################

par(mar=c(4, 2.5, 3, 1), mgp = c(2.5, 1, 0))
par(las=1,lwd=1,cex.main=font, cex.lab=1.6, cex.axis=1.3)
par(cex.main=1.6)
par(las=1)
zehner=0
for (i in 1:9){
  plot(output$relfreq.hist[10*zehner+i,],pch=20,xlab="Iteration",ylab="",main=paste("Proposal probability for X",10*zehner+i, sep=""),ylim=c(0,1),lwd=2, type = "l")
  #points(output_PO$relfreq.hist[10*zehner+i,],pch=20,xlab="Iteration",ylab="",main=paste("Var",10*zehner+i),ylim=c(0,1),col="yellow",cex=0.3)
  lines(output2$relfreq.hist[10*zehner+i,],pch=20,xlab="Iteration",ylab="",main=paste("Proposal probability for X",10*zehner+i, sep=""),ylim=c(0,1),col="orange",lwd=2)
  lines(output3$relfreq.hist[10*zehner+i,],pch=20,xlab="Iteration",ylab="",main=paste("Proposal probability for X",10*zehner+i, sep=""),ylim=c(0,1),col="purple",lwd=2)
  abline(h=out$probne0[-1][i],col="red", lwd =2, lty=2) 
}

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend = c("L = p", "L = p/n", "L = 100p", "True posterior"), 
       col=c("black", "orange", "purple", "red"), cex=2, box.lty=1, 
      pt.bg = 'white', lty=c(1,1,1,2), lwd =2, ncol=4) #,
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

