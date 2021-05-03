par(las=1)

load("./Files/R_simulation_saves/MAdaSub_results_Golub_random_q_2_5_L_05_2_p_epsilon.RData")


# computational time
results$time

# acceptance rate
summary(colMeans(results$acc[21:50,1:50] / 50000))

results$time 
results$acc[,2] 

repeats = dim(results$importance.individual)[1] - 1
nb.k = dim(results$importance.individual)[2]
frac.par.update = 0.5

parameter.list = results$parameter.list 

importance.individual = results$importance.individual 

indices = c()
for (i in 1:nb.k)
  indices = sort(unique(c(indices, which(importance.individual[repeats+1,i,]>0.1))))


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
Xnames = golub_names 
#############################################################

indices
golub_names[indices]
summary(importance.individual[repeats+1,1:25,956])
summary(importance.individual[repeats+1,1:25,2481])


#indices = which(importance.individual[repeats+1,50,]>0.05)

# golub_names[indices]

par(mfrow = c(2,4))
par(cex.main=1.5)
size_axis <- 1
par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))

boxplot(importance.individual[1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 1 - 50,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[2,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 50,001 - 100,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[3,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 100,001 - 150,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[repeats+1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 1,000,001 - 2,500,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)

boxplot(importance.individual[1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 1 - 50,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[2,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 50,001 - 100,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[3,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 100,001 - 150,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[repeats+1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 1,000,001 - 2,500,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)


########################################################################
## Evolution of empirical inclusion probabilities 

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
  plot(y_median, ylim=c(0,1), xlim=c(0,nb.k),
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

mtext("*each round consists of 50,000 iterations", side=1, adj=1, line = 4.75 )




