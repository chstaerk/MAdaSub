library("BAS")
library("coda")

#setwd("C:\\Users\\Staerk\\ScieboBonn\\MAdaSub PAPER\\R files")
#source("MAdaSub_LOAD.R")
RNGkind(kind = NULL, normal.kind = NULL, sample.kind = NULL)

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

##################################################################
# Begin MAdaSub algo
##################################################################

family="normal"

savings = 1

q=5

epsilon = 1/p # 1e-7
priormean = q/p # rep(q/p,p)
L = p #rep(p,p)  #K=p/L
#L = 1000
#priorprob= 1/2
#g=n
priorprob= q/p

#prior="g" # g-prior (original form with Jeffreys prior on variance)
#prior="EBIC"
prior = "independence"

const=0.6
g = 5 # for Trecator data 
hyper = FALSE
a_prior = 1  
b_prior = (p-5)/5 # (p-5)/5

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)


choices <-c("290,000 iterations (as in Supplement)", "Custom number of iterations (should be a multiple of 10)")
mychoice <- menu( choices , graphics=TRUE, title="How many iterations per chain?" )

if (mychoice==1) {
  Iter = c(rep(5000,20),190000)
  n_Iter = 290000
}
if (mychoice==2) {
  n_Iter = -1
  while(n_Iter < 1 ){
    n_Iter <- readline("Enter the number of iterations: ")
    n_Iter <- ifelse(grepl("\\D",n_Iter),-1,as.integer(n_Iter))
    if(is.na(n_Iter)){break}  # breaks when hit enter
  }
  Iter = c(rep(n_Iter/10, 5), n_Iter/2)
}


#Iter = c(rep(5000,20),190000)
set.seed(22)
start.time <- Sys.time()
output=MAdaSub_wrapper(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,
                       plot=FALSE)
end.time <- Sys.time()
end.time - start.time
comp_time = as.numeric(difftime(end.time, start.time, units = "secs"))

###
run_2Mio =FALSE # set to TRUE if MAdaSub algorithm should also be run for 2,000,000 iterations 

if (run_2Mio) {
  Iter = c(100000,1900000)
  set.seed(22)
  start.time <- Sys.time()
  output_2Mio=MAdaSub_wrapper(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,
                       plot=FALSE)
  end.time <- Sys.time()
  end.time - start.time
  comp_time_2Mio = as.numeric(difftime(end.time, start.time, units = "secs"))
}

#setwd("C://Users//Staerk//sciebo//MAdaSub PAPER//R Simulation Saves")
#save(output, file="MAdaSub_Highdim_Tecator_Lamnisos_190K_with100Kburnin_q5.RData")

#########################################################################
#load("MAdaSub_Highdim_Golub_400K_q10.RData") # output


dev.off()
win.graph(width = 10, height = 6, pointsize = 10)
plot(output$relfreq.final,pch=20,ylab="Final proposal probabilities", xlab="Covariate index", main=paste("MAdaSub proposal probabilities for Tecator data after", n_Iter, "iterations"), ylim=c(0,1))


##### Compute effective sample size #############################
Gamma_indi_mat = output$CS
for (i in (2:(dim(output$CS)[2]))) {
  Gamma_indi_mat[,i] = output$CS[,i] - output$CS[,i-1]
}
Gamma_indi_mcnc= as.mcmc(t(Gamma_indi_mat)) 
eff_size = effectiveSize(Gamma_indi_mcnc)
effective_size_median = median(eff_size)

effective_size_median 
print(paste("Median effective sample size (for last", Iter[length(Iter)],"iterations):", round(effective_size_median)))

# compute effective sample sizes for 2,000,000 iterations with thinning every 10th iterations and 100,000 burn-in iterations
if (run_2Mio) {  
  indices_thinned <- seq(10,1900000,10)
  Gamma_indi_mat = matrix(0, 100, length(indices_thinned) )
  for (i in (2:length(indices_thinned))) {
    Gamma_indi_mat[,i] = output_2Mio$CS[,indices_thinned[i]] - output_2Mio$CS[,indices_thinned[i]-1]
  }
  Gamma_indi_mcnc= as.mcmc(t(Gamma_indi_mat)) 
  eff_size = effectiveSize(Gamma_indi_mcnc)
  effective_size_median_thinning = median(eff_size)

  effective_size_median_thinning
}


##################################################################################
### Additional plots 
##################################################################################

#font = 1
#par(las=1,lwd=1,cex.main=font, cex.lab=font, cex.axis=font)
#par(cex.main=1.5)
#options(scipen=5)

#plot(output$values,pch=20,xlab="Iteration",ylab="",main="Evolution of sampled values")

#par(mfrow=c(1,1))
#plot(output$acc.prob,pch=20,ylim=c(0,1),xlab="Iteration",ylab="", main="Acceptance rate")


### Plot of sampled sizes of V^(t) and f(V^(t))
#par(mfrow=c(1,1))
#plot(output$V.size,pch=20,ylim=c(0,max(output$V.size)),xlab="Iteration",ylab="Model size",main="Evolution of model sizes")
#points(output$S.size,col='red',pch=20,cex=0.7)












