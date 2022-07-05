library("MASS")
#library("Matrix")
library("RcppNumerical")
#library("BMS")
#library("BAS")

library("parallel")

#source("MAdaSub_LOAD.R")

######### Golub data ##################################
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

family = "binomial"

const=1
savings = 100

q=2
epsilon = 1/p # 1/(2*p) # 1e-7
priormean = rep(q/p,p)
L = p #rep(p,p)  #K=p/L
priorprob = q/p # 1/2

prior="EBIC"
#prior = "independence"

#g = 5 # for Trecator data 
g = 100 # for PCR and Golub data 
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
#numCores = 50 

 
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

#nb.k = numCores

 choices <-c( "50 rounds (as in paper)                                          ", 
              "Custom number of rounds                                          ")
 mychoice <- menu( choices , graphics=TRUE, title="How many rounds for each chain (including burn-in)?" )
 
 if (mychoice==1) {
   repeats=50 
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

 choices <-c( "20,000 iterations per round (as in paper)                      ", 
              "Custom number of iterations per round                        ")
 mychoice <- menu( choices , graphics=TRUE, title="How many iterations per round for each chain?" )
 
 if (mychoice==1) {
   M = 20000
 }
 if (mychoice==2) {
   M = -1
   while(M < 1 ){
     M <- readline("Enter the number of iterations per round: ")
     M <- ifelse(grepl("\\D",M),-1,as.integer(M))
     if(is.na(M)){break}  # breaks when hit enter
   }
 }

frac.par.update = 0.5 # fraction of chains with parallel updating
# if frac.par.update = 1, all chains used for parallel updating
# if frac.par.update = 0, independent chains (no parallel updating)

###
#numCores = 1
#nb.k = 2
###
RNGkind("L'Ecuyer-CMRG")
set.seed(21)

#nb.k = 50 # 50
#numCores = 5 # 50 # nb.k Change to 50 for comp time !
#repeats = 50 # 50
#burnin_rounds = 20 # 20
#M = 20000 # 50000

###
priormean <- list()
for (i in 1:nb.k) {
  priormean[[i]] <- rep(runif(1,2,5)/p, p)
}

L <- list()
for (i in 1:nb.k) {
  L[[i]] <- runif(1, p/2, 2*p)
}
###

start.time <- Sys.time()
results = MAdaSub_parallel(data=data,nb.k=nb.k,M=M,repeats=repeats,priormean=priormean,L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=FALSE,a_prior=NULL,b_prior=NULL,
                           frac.par.update=frac.par.update, burnin_rounds = burnin_rounds, numCores = numCores) 
end.time <- Sys.time()
time = as.double(end.time-start.time, units = "secs")

results$time = time
save(results, file="MAdaSub_results_Golub_random_q_2_5_L_05_2_p_epsilon20K.RData")


#results = mclapply(parameter.list, f_parallel, mc.cores = numCores, 
#                   data=data ,M=M ,const=const ,savings=savings ,family=family ,epsilon=epsilon ,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,
#                   mc.set.seed = TRUE,  mc.preschedule = FALSE )

 
 parameter.list = results$parameter.list 
 
 importance.individual = results$importance.individual 
 indices = c()
 for (i in 1:nb.k)
   indices = sort(unique(c(indices, which(importance.individual[repeats+1,i,]>0.1))))
 
 
 dev.off()
 win.graph(width = 14, height = 7, pointsize = 10)
 
 size_axis <- 0.9
 
 par(mfrow = c(2,4))
 par(cex.main=1.5)
 par(las=1)
 par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))
 
 boxplot(importance.individual[1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
         main=paste("Serial updating \n Iterations ", 1 ," - ", M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
 axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
 axis(1, at = 1:100, labels=rep("",100))
 boxplot(importance.individual[2,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
         main=paste("Serial updating \n Iterations ", M+1 ," - ", 2*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
 axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
 axis(1, at = 1:100, labels=rep("",100))
 boxplot(importance.individual[3,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
         main=paste("Serial updating \n Iterations ", 2*M+1 ," - ", 3*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
 axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
 axis(1, at = 1:100, labels=rep("",100))
 boxplot(importance.individual[repeats+1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
         main=paste("Serial updating \n Iterations ", burnin_rounds*M+1 ," - ", repeats*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
 axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
 axis(1, at = 1:100, labels=rep("",100))
 #points(results$parameter.list[[1]]$priormean.cur, col="red",pch=20)
 
 boxplot(importance.individual[1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
         main=paste("Parallel updating \n Iterations ", 1 ," - ", M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
 axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
 axis(1, at = 1:100, labels=rep("",100))
 boxplot(importance.individual[2,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
         main=paste("Parallel updating \n Iterations ", M+1 ," - ", 2*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
 axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
 axis(1, at = 1:100, labels=rep("",100))
 boxplot(importance.individual[3,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
         main=paste("Parallel updating \n Iterations ", 2*M+1 ," - ", 3*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
 axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
 axis(1, at = 1:100, labels=rep("",100))
 boxplot(importance.individual[repeats+1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
         main=paste("Parallel updating \n Iterations ", burnin_rounds*M+1 ," - ", repeats*M, sep=""), ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
 axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
 axis(1, at = 1:100, labels=rep("",100))
 #points(results$parameter.list[[1]]$priormean.cur, col="red",pch=20)
 ###################################################################################################
 
 Xnames = golub_names 
 indices = sort(order(importance.individual[repeats+1,1,], decreasing=TRUE)[1:4])
 
 win.graph(width = 14, height = 7, pointsize = 10)
 #par(mfrow = c(2,4))
 par(cex.main=1.2)
 
 ############################################################################
 m <- matrix(c(1,2,3,4,5,6,7,8,9,9,9,9),nrow = 3,ncol = 4,byrow = TRUE)
 layout(mat = m, heights = c(0.48,0.48,0.04))
 ############################################################################
 
 size_axis <- 1
 
 #ayout(matrix(c(1:16), nrow = 2, ncol = 8, byrow = TRUE), widths=rep(c(1/4*4/5, 1/4*1/5), 4), heights=c(0.5, 0.5))
 
 for (i in 1:4) {
   par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))
   y_median <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,median)
   y_q05 <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,quantile,probs=0.05)
   y_q95 <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,quantile,probs=0.95)
   y_min <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,min)
   y_max <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,max)
   plot(y_median, ylim=c(0,1), 
        main=paste("Serial updating \n Variable ", indices[i], " (gene ", Xnames[indices[i]], ")", sep=""), ylab="Emp. inclusion freq.", xlab="Number of rounds*", cex.axis = size_axis,type="l",cex.lab=1.2, cex.axis=1.05)
   lines(y_q05)
   lines(y_q95)
   polygon(c(1:repeats,rev(1:repeats)),c(y_q05,rev(y_q95)),col="gray")
   #lines(y_min)
   #lines(y_max)
   #polygon(c(1:nb.k,rev(1:nb.k)),c(y_min,rev(y_max)),col="gray")
   lines(y_median, lwd=2)
 }
 
 for (i in 1:4) {
   par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))
   y_median <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,median)
   y_q05 <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,quantile,probs=0.05)
   y_q95 <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,quantile,probs=0.95)
   y_min <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,min)
   y_max <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,max)
   plot(y_median, ylim=c(0,1),
        main=paste("Parallel updating \n Variable ", indices[i], " (gene ", Xnames[indices[i]], ")", sep=""), ylab="Emp. inclusion freq.", xlab="Number of rounds*", cex.axis = size_axis,type="l",cex.lab=1.2, cex.axis=1.05)
   lines(y_q05)
   lines(y_q95)
   polygon(c(1:repeats,rev(1:repeats)),c(y_q05,rev(y_q95)),col="gray")
   #lines(y_min)
   #lines(y_max)
   #polygon(c(1:nb.k,rev(1:nb.k)),c(y_min,rev(y_max)),col="gray")
   lines(y_median, lwd=2)
 }
 
 mtext(paste("*each round consists of ", M," iterations", sep=""), side=1, adj=1, line = 4.75 )
 
 
