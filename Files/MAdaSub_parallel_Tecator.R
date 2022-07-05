library("MASS")
library("parallel")


#source("MAdaSub_LOAD.R")

######### Tecator data ##################################
library("fda.usc")
data(tecator)
data = list()
samples = 1:172
data$x = tecator$absorp.fdata$data[samples,]
data$x = scale(data$x,scale = FALSE)
#data$x = t(scale(t(data$x)))
data$y = tecator$y[["Fat"]][samples]
n = length(samples)
p = dim(data$x)[2]
family="normal"
#############################################################


const=1
savings = 100

q=5
epsilon = 1/p # 1/(2*p) # 1/(2*p) # 1e-7
priormean = rep(q/p,p)
L = p #rep(p,p)  #K=p/L
priorprob = q/p # 1/2
#priorprob = 1/2

#prior="EBIC"
prior = "independence"

g = 5 # for Trecator data 
#g = 100 # for PCR data 
#g = p^2 # p^2
hyper = FALSE
a_prior = 1  
b_prior = (p-5)/5 # (p-5)/5



choices <-c( "50 CPUs (as in paper)                        ", 
             "Custom number of CPUs (for Windows use 1 CPU)")
mychoice <- menu( choices , graphics=TRUE, title="How many CPUs?" )

if (mychoice==1) {
  numCores=50 
}
if (mychoice==2) {
  numCores = -1
  while(numCores < 1 ){
    numCores <- readline("Enter the number of CPUs (for Windows use 1 CPU): ")
    numCores <- ifelse(grepl("\\D",numCores),-1,as.integer(numCores))
    if(is.na(numCores)){break}  # breaks when hit enter
  }
}

#numCores <- detectCores()

choices <-c( "50 chains (as in paper)                    ", 
             "Custom number of chains                    ")
mychoice <- menu( choices , graphics=TRUE, title="How many MAdaSub chains?" )

if (mychoice==1) {
  nb.k=50 
}
if (mychoice==2) {
  nb.k = -1
  while(nb.k < 1 ){
    nb.k <- readline("Enter the number of chains (minimum 4 chains): ")
    nb.k <- ifelse(grepl("\\D",nb.k),-1,as.integer(nb.k))
    if(is.na(nb.k)){break}  # breaks when hit enter
  }
}

choices <-c( "58 rounds (as in paper)                                          ", 
             "Custom number of rounds                                          ")
mychoice <- menu( choices , graphics=TRUE, title="How many rounds for each chain (including burn-in)?" )

if (mychoice==1) {
  repeats=58 
}
if (mychoice==2) {
  repeats = -1
  while(repeats < 1 ){
    repeats <- readline("Enter the number of rounds per chain (minimum 4 rounds): ")
    repeats <- ifelse(grepl("\\D",repeats),-1,as.integer(repeats))
    if(is.na(repeats)){break}  # breaks when hit enter
  }
}

choices <-c( "20 burn-in rounds (as in paper)                           ", 
             "Custom number of burn-in rounds                           ")
mychoice <- menu( choices , graphics=TRUE, title="How many burn-in rounds for each chain?" )

if (mychoice==1) {
  burnin_rounds=20 
}
if (mychoice==2) {
  burnin_rounds = -1
  while(burnin_rounds < 1 ){
    burnin_rounds <- readline(paste("Enter the number of burn-in rounds per chain (between", 1 ,"and", repeats-2, "rounds):"))
    burnin_rounds <- ifelse(grepl("\\D",burnin_rounds),-1,as.integer(burnin_rounds))
    if(is.na(burnin_rounds)){break}  # breaks when hit enter
  }
}


choices <-c( "5000 iterations per round (as in paper)                      ", 
             "Custom number of iterations per round                        ")
mychoice <- menu( choices , graphics=TRUE, title="How many iterations per round for each chain?" )

if (mychoice==1) {
  M = 5000
}
if (mychoice==2) {
  M = -1
  while(M < 1 ){
    M <- readline("Enter the number of iterations per round: ")
    M <- ifelse(grepl("\\D",M),-1,as.integer(M))
    if(is.na(M)){break}  # breaks when hit enter
  }
}

#nb.k = 50
#numCores = 50 
#repeats = 58
#M = 5000 #  
#burnin_rounds = 20


frac.par.update = 0.5 # fraction of chains with parallel updating
# if frac.par.update = 1, all chains used for parallel updating
# if frac.par.update = 0, independent chains (no parallel updating)

RNGkind("L'Ecuyer-CMRG")
set.seed(21)

###
priormean <- list()
for (i in 1:nb.k) {
  priormean[[i]] <- rep(runif(1,2,10)/p, p)
}

L <- list()
for (i in 1:nb.k) {
  L[[i]] <- runif(1, p/2, 2*p)
}
#L = p
###

set.seed(1)

start.time <- Sys.time()
results = MAdaSub_parallel(data=data,nb.k=nb.k,M=M,repeats=repeats,priormean=priormean,L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=FALSE,a_prior=NULL,b_prior=NULL,
                           frac.par.update=frac.par.update, burnin_rounds=burnin_rounds, numCores=numCores) 
end.time <- Sys.time()
time = as.double(end.time-start.time, units = "secs")

results$time = time
save(results, file="MAdaSub_results_Tecator_random_q_2_10_L_05_2p_epsilon.RData")


parameter.list = results$parameter.list 

importance.individual = results$importance.individual 
indices = 1:100

dev.off()
win.graph(width = 14, height = 7, pointsize = 10)

size_axis <- 0.9

par(mfrow = c(2,4))
par(cex.main=1.5)
par(las=1)
par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))

boxplot(importance.individual[1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main=paste("Serial updating \n Iterations ", 1 ," - ", M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[2,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main=paste("Serial updating \n Iterations ", M+1 ," - ", 2*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[3,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main=paste("Serial updating \n Iterations ", 2*M+1 ," - ", 3*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[repeats+1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main=paste("Serial updating \n Iterations ", burnin_rounds*M+1 ," - ", repeats*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
#points(results$parameter.list[[1]]$priormean.cur, col="red",pch=20)

boxplot(importance.individual[1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main=paste("Parallel updating \n Iterations ", 1 ," - ", M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[2,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main=paste("Parallel updating \n Iterations ", M+1 ," - ", 2*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[3,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main=paste("Parallel updating \n Iterations ", 2*M+1 ," - ", 3*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[repeats+1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main=paste("Parallel updating \n Iterations ", burnin_rounds*M+1 ," - ", repeats*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
#points(results$parameter.list[[1]]$priormean.cur, col="red",pch=20)
###################################################################################################


