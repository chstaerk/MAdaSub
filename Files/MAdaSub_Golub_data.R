library("BAS")
library("coda")

#setwd("C:\\Users\\Staerk\\ScieboBonn\\MAdaSub PAPER\\R files")
#source("MAdaSub_LOAD.R")

######### Golub DATA ##################################
data=list()

library("golubEsets")
library("genefilter")

data("Golub_Train")
data("Golub_Merge")

phenoTrain <- phenoData(Golub_Merge)
n= dim(exprs(Golub_Merge))[2]
Y = phenoTrain[[2]]
data$y=numeric(n)
data$y[Y=="ALL"] = 1  # 1 == ALL, 0 == AML
X <- exprs(Golub_Merge)

#Thresholding
dim(X)
X[X < 100] <- 100
X[X > 16000] <- 16000

#Filtering:
mmfilt <- function(r = 5, d = 500, na.rm = TRUE) {
  function(x) {
    minval <- min(x, na.rm = na.rm)
    maxval <- max(x, na.rm = na.rm)
    (maxval/minval > r) && (maxval - minval > d)
  }
}
mmfun <- mmfilt()
ffun <- filterfun(mmfun)
good <- genefilter(X, ffun)
sum(good)
X <- X[good, ]

#Log-transform
X <- log10(X)

#Normalize 
X <- scale(X)

data$x = t(X)
data$x = scale(data$x, scale=FALSE) # mean-center the covariates
p=dim(data$x)[2]
golub_names = colnames(data$x)
#############################################################


##################################################################
# Begin MAdaSub algo
##################################################################
family="binomial"  

savings = 100

### Tuning parameters for MAdaSub ####
q = 5
priormean = q/p # prior expected search size (for MAdaSub algorithm) ### rep(q/p,p)
L = p #rep(p,p)  # corresponding to K=p/L in AdaSub 
epsilon = 1/p 

### Type of prior ####
prior="EBIC" # posterior resulting from EBIC approximation 

### Prior parameters #####
const = 1 # constant for EBIC 

choices <-c("1,000,000 iterations (as in paper)", "Custom number of iterations (should be a multiple of 50)")
mychoice <- menu( choices , graphics=TRUE, title="How many iterations per chain?" )

if (mychoice==1) {
  Iter = c(rep(20000,50)) 
}
if (mychoice==2) {
  n_Iter = -1
  while(n_Iter < 1 ){
    n_Iter <- readline("Enter the number of iterations: ")
    n_Iter <- ifelse(grepl("\\D",n_Iter),-1,as.integer(n_Iter))
    if(is.na(n_Iter)){break}  # breaks when hit enter
  }
  Iter = c(rep(n_Iter/50, 50))
}


#Iter = c(rep(50000,50)) # vector of number of iterations

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)

q = 10
priormean = q/p 

print("First chain with q=10")
set.seed(1)
start.time <- Sys.time()
output=MAdaSub_wrapper(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior, plot=FALSE)
#output=MAdaSub(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior, plot=FALSE)
end.time <- Sys.time()
end.time - start.time
comp_time = as.numeric(difftime(end.time, start.time, units = "secs"))

q=5
priormean = q/p 

print("Second chain with q=5")
set.seed(2)
start.time <- Sys.time()
output2=MAdaSub_wrapper(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior, plot=FALSE)
#output2=MAdaSub(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior, plot=FALSE)
end.time <- Sys.time()
end.time - start.time
comp_time2 = as.numeric(difftime(end.time, start.time, units = "secs"))

q=2
priormean = q/p 

print("Third chain with q=2")
set.seed(3)
start.time <- Sys.time()
output3=MAdaSub_wrapper(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior, plot=FALSE)
#output3=MAdaSub(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior, plot=FALSE)
end.time <- Sys.time()
end.time - start.time
comp_time3 = as.numeric(difftime(end.time, start.time, units = "secs"))


#getwd()
#setwd("C://Users//Staerk//scieboBonn//MAdaSub PAPER//R Simulation Saves")
save(output, file="MAdaSub_Highdim_Golub_1Mio_q10.RData")
save(output2, file="MAdaSub_Highdim_Golub_1Mio_q5.RData")
save(output3, file="MAdaSub_Highdim_Golub_1Mio_q2.RData")


#########################################################################
#load("./Files/R_simulation_saves/MAdaSub_Highdim_Golub_1Mio_q10.RData") # output
#load("./Files/R_simulation_saves/MAdaSub_Highdim_Golub_1Mio_q5.RData")  # output2
#load("./Files/R_simulation_saves/MAdaSub_Highdim_Golub_1Mio_q2.RData")  # output3

#which(output$relfreq.final>=0.5) # Median Probability model
#output$best.S

############################
#golub_names[which(output$relfreq.final>=0.2)]
#output$relfreq.final[output$relfreq.final>=0.2]
#output2$relfreq.final[output$relfreq.final>=0.2]
#output3$relfreq.final[output$relfreq.final>=0.2]
############################


#############################################
## Plots for paper
#############################################
dev.off()
win.graph(width = 12, height = 6, pointsize = 10)
par(cex.main=1.5)
options(scipen=5)
par(las=1,lwd=1, cex.lab=1.6, cex.axis=1.3)
par(cex.main=1.6)

# Scatterplots: Estimated marginal incl prob from different runs
par(mfrow=c(1,3))
plot(output3$relfreq.final,output2$relfreq.final,pch=20,xlab="Run 1",ylab="Run 2",xlim=c(0,1),ylim=c(0,1),cex=1.5,main="q = 2 vs. q = 5")
abline(a=0,b=1)
plot(output3$relfreq.final,output$relfreq.final,pch=20,xlab="Run 1",ylab="Run 3",xlim=c(0,1),ylim=c(0,1),cex=1.5,main="q = 2 vs. q = 10")
abline(a=0,b=1)
plot(output2$relfreq.final,output$relfreq.final,pch=20,xlab="Run 2",ylab="Run 3",xlim=c(0,1),ylim=c(0,1),cex=1.5,main="q = 5 vs. q = 10")
abline(a=0,b=1)


#############################################################################################
##### MAdaSub on Leukemia data with only 200,000 iterations for plots (same seeds as before) #####
#############################################################################################

run_200K <- FALSE

if (run_200K) {

family="binomial"  

savings = 100 

### Tuning parameters for MAdaSub ####
q = 10
priormean = q/p # prior expected search size (for MAdaSub algorithm) ### rep(q/p,p)
L = p #rep(p,p)  # corresponding to K=p/L in AdaSub 
epsilon = 1/p #1e-7

### Type of prior ####
prior="EBIC" # posterior resulting from EBIC approximation 

### Prior parameters #####
const = 1 # constant for EBIC 

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)


Iter = c(200000) # vector of number of iterations

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)

q = 10
priormean = q/p 

set.seed(1)
start.time <- Sys.time()
output=MAdaSub(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,plot=FALSE)
end.time <- Sys.time()
end.time - start.time
comp_time = as.numeric(difftime(end.time, start.time, units = "secs"))

q=5
priormean = q/p 

set.seed(2)
start.time <- Sys.time()
output2=MAdaSub(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,plot=FALSE)
end.time <- Sys.time()
end.time - start.time
comp_time2 = as.numeric(difftime(end.time, start.time, units = "secs"))

q=2
priormean = q/p 

set.seed(3)
start.time <- Sys.time()
output3=MAdaSub(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,plot=FALSE)
end.time <- Sys.time()
end.time - start.time
comp_time3 = as.numeric(difftime(end.time, start.time, units = "secs"))

par(mfrow=c(1,3))
plot(output$relfreq.final,output2$relfreq.final,pch=20,xlab="Run 1",ylab="Run 2",xlim=c(0,1),ylim=c(0,1),cex=1.5)
abline(a=0,b=1)
plot(output$relfreq.final,output3$relfreq.final,pch=20,xlab="Run 1",ylab="Run 3",xlim=c(0,1),ylim=c(0,1),cex=1.5)
abline(a=0,b=1)
plot(output2$relfreq.final,output3$relfreq.final,pch=20,xlab="Run 2",ylab="Run 3",xlim=c(0,1),ylim=c(0,1),cex=1.5)
abline(a=0,b=1)


#setwd("C://Users//Staerk//scieboBonn//MAdaSub PAPER//R Simulation Saves")
#save(output, file="MAdaSub_Highdim_Golub_200new_q10.RData")
#save(output2, file="MAdaSub_Highdim_Golub_200new_q5.RData")
#save(output3, file="MAdaSub_Highdim_Golub_200new_q2.RData")

#########################################################################
#load("./Files/R_simulation_saves/MAdaSub_Highdim_Golub_200new_q10.RData") # output
#load("./Files/R_simulation_saves/MAdaSub_Highdim_Golub_200new_q5.RData")  # output2
#load("./Files/R_simulation_saves/MAdaSub_Highdim_Golub_200new_q2.RData")  # output3


#############################################
## Plots 
#############################################

par(cex.main=1.5)
options(scipen=5)
par(las=1,lwd=1, cex.lab=1.6, cex.axis=1.3)
par(cex.main=1.6)


### Plot of sampled sizes of V^(t) and f(V^(t))
par(mfrow=c(3,1))
plot(output3$V.size,pch=20,ylim=c(0,20),xlab="Iteration",ylab="Model size",main="q = 2")
points(output3$S.size,col='red',pch=20,cex=0.7)
plot(output2$V.size,pch=20,ylim=c(0,20),xlab="Iteration",ylab="Model size",main="q = 5")
points(output2$S.size,col='red',pch=20,cex=0.7)
plot(output$V.size,pch=20,ylim=c(0,max(output$V.size)),xlab="Iteration",ylab="Model size", main="q = 10")
points(output$S.size,col='red',pch=20,cex=0.7)

}



