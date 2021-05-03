library("MASS")
library("parallel")


source("MAdaSub_LOAD.R")

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



numCores <- detectCores()
numCores = 50 #70
nb.k = numCores
repeats = 58
M = 5000 #  for Tecator data: M = 5000, for PCR data: M = 20000
burnin_rounds = 20

frac.par.update = 0.5 # fraction of chains with parallel updating
# if frac.par.update = 1, all chains used for parallel updating
# if frac.par.update = 0, independent chains (no parallel updating)

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

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

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

#start.time <- Sys.time()
#results = MAdaSub_parallel(data=data,nb.k=nb.k,M=M,repeats=repeats,priormean=priormean,L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=FALSE,a_prior=NULL,b_prior=NULL,
#                           frac.par.update=frac.par.update, burnin_rounds=burnin_rounds) 
#end.time <- Sys.time()
#time = as.double(end.time-start.time, units = "secs")

priormean <- priormean[[1]]
L <- L[[1]]

Iter = c(rep(5000,58))
set.seed(22)
start.time <- Sys.time()
output=MAdaSub_wrapper(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,
                       plot=FALSE)
end.time <- Sys.time()
end.time - start.time
comp_time = as.numeric(difftime(end.time, start.time, units = "secs"))

save(comp_time, file="MAdaSub_serial_timing_Tector.RData")

#results = mclapply(parameter.list, f_parallel, mc.cores = numCores, 
#                   data=data ,M=M ,const=const ,savings=savings ,family=family ,epsilon=epsilon ,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,
#                   mc.set.seed = TRUE,  mc.preschedule = FALSE )


