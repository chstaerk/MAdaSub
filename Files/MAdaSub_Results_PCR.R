par(las=1)

load("./Files/R_simulation_saves/MAdaSub_results_PCR_random_q_2_5_L_05_2_p_epsilon20K.RData")


# computational time
results$time

# acceptance rate
summary(colMeans(results$acc[21:50,1:50] / 20000))
#summary(colMeans(results$acc[21:50,1:50] / 50000))
#summary(colMeans(results$acc[21:50,1:25] / 50000))
#summary(colMeans(results$acc[21:50,26:50] / 50000))

repeats = dim(results$importance.individual)[1] - 1
nb.k = dim(results$importance.individual)[2]
frac.par.update = 0.5

parameter.list = results$parameter.list 

importance.individual = results$importance.individual 

indices = c()
for (i in 1:nb.k)
  indices = sort(unique(c(indices, which(importance.individual[repeats+1,i,]>0.05))))

Xnames = scan("Files/gene_id.txt",what=character())
indices
Xnames[indices]
summary(importance.individual[repeats+1,1:50,7639])
summary(importance.individual[repeats+1,1:50,7640])


par(mfrow = c(2,4))
par(cex.main=1.5)

size_axis <- 1
par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))

boxplot(importance.individual[1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 1 - 20,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[2,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 20,001 - 40,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[3,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 40,001 - 60,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[repeats+1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 200,001 - 1,000,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)

boxplot(importance.individual[1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 1 - 20,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[2,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 20,001 - 40,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[3,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 40,001 - 60,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)
boxplot(importance.individual[repeats+1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 200,001 - 1,000,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis)



########################################################################
## Evolution of empirical inclusion probabilities 

#par(mfrow = c(2,5))
par(cex.main=1.31)
par(cex.lab=1.45)

size_axis <- 1.2

############################################################################
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,11,11,11,11),nrow = 3,ncol = 5,byrow = TRUE)
layout(mat = m, heights = c(0.48,0.48,0.04))
############################################################################

for (i in 1:5) {
  par(mar=c(4, 4, 4, 1), mgp = c(2.7, 1, 0))
  y_median <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,median)
  y_q05 <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,quantile,probs=0.05)
  y_q95 <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,quantile,probs=0.95)
  y_min <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,min)
  y_max <- apply(importance.individual[1:repeats,(frac.par.update*nb.k+1):nb.k,indices[i]],1,max)
  plot(y_median, ylim=c(0,1), xlim=c(0,nb.k),
       main=paste("Serial updating \n Variable ", indices[i], " (gene ", Xnames[indices[i]], ")", sep=""), ylab="Emp. inclusion freq.", xlab="Number of rounds*", cex.axis = size_axis,type="l")
  lines(y_q05)
  lines(y_q95)
  polygon(c(1:repeats,rev(1:repeats)),c(y_q05,rev(y_q95)),col="gray")
  #lines(y_min)
  #lines(y_max)
  #polygon(c(1:nb.k,rev(1:nb.k)),c(y_min,rev(y_max)),col="gray")
  lines(y_median, lwd=2)
}

for (i in 1:5) {
  par(mar=c(4, 4, 4, 1), mgp = c(2.7, 1, 0))
  y_median <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,median)
  y_q05 <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,quantile,probs=0.05)
  y_q95 <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,quantile,probs=0.95)
  y_min <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,min)
  y_max <- apply(importance.individual[1:repeats,1:(frac.par.update*nb.k),indices[i]],1,max)
  plot(y_median, ylim=c(0,1), xlim=c(0,nb.k),
       main=paste("Parallel updating \n Variable ", indices[i], " (gene ", Xnames[indices[i]], ")", sep=""), ylab="Emp. inclusion freq.", xlab="Number of rounds*", cex.axis = size_axis,type="l")
  lines(y_q05)
  lines(y_q95)
  polygon(c(1:repeats,rev(1:repeats)),c(y_q05,rev(y_q95)),col="gray")
  #lines(y_min)
  #lines(y_max)
  #polygon(c(1:nb.k,rev(1:nb.k)),c(y_min,rev(y_max)),col="gray")
  lines(y_median, lwd=2)
}

mtext("*each round consists of 20,000 iterations", side=1, adj=1, line = 5, cex=1.2)





###########################################################







