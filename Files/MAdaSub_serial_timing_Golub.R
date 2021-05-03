library("MASS")
#library("Matrix")
library("RcppNumerical")
#library("BMS")
#library("BAS")

library("parallel")

source("MAdaSub_LOAD.R")

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



#numCores <- detectCores()
numCores = 50 #70
nb.k = numCores
repeats = 50 # 500 # 40
M = 50000 # 20000 #  for Tecator data: M = 5000, for Golub and PCR data: M = 20000
burnin_rounds = 20

frac.par.update = 0.5 # 0.5 # fraction of chains with parallel updating
# if frac.par.update = 1, all chains used for parallel updating
# if frac.par.update = 0, independent chains (no parallel updating)

###
#numCores = 1
#nb.k = 2
###
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

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

priormean <- priormean[[1]]
L <- L[[1]]

Iter = c(rep(50000,50))
set.seed(22)
start.time <- Sys.time()
output=MAdaSub_wrapper(data,Iter=Iter,priormean=priormean, L=L,const=const,savings=savings,family=family,epsilon=epsilon,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,
                       plot=FALSE)
end.time <- Sys.time()
end.time - start.time
comp_time = as.numeric(difftime(end.time, start.time, units = "secs"))

save(comp_time, file="MAdaSub_serial_timing_Golub.RData")



#results = mclapply(parameter.list, f_parallel, mc.cores = numCores, 
#                   data=data ,M=M ,const=const ,savings=savings ,family=family ,epsilon=epsilon ,priorprob=priorprob,prior=prior,g=g,hyper=hyper,a_prior=a_prior,b_prior=b_prior,
#                   mc.set.seed = TRUE,  mc.preschedule = FALSE )


