par(las=1)

load("./Files/R_simulation_saves/MAdaSub_results_Tecator_random_q_2_10_L_05_2p_epsilon.RData")

M = 5000

results$time

repeats = dim(results$importance.individual)[1] - 1
nb.k = dim(results$importance.individual)[2]
frac.par.update = 0.5

# mean acceptance rate after burnin period of 100,000 iterations
mean(results$acc[repeats,]/5000)


parameter.list = results$parameter.list 

importance.individual = results$importance.individual 
indices = 1:100

size_axis <- 0.9

par(mfrow = c(2,4))
par(cex.main=1.5)
par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))

boxplot(importance.individual[1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 1 - 5,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[2,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 5,001 - 10,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[3,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 10,001 - 15,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[repeats+1,(frac.par.update*nb.k+1):nb.k,indices], ylim=c(0,1), names=indices,
        main="Serial updating \n Iterations 100,001 - 290,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
#points(results$parameter.list[[1]]$priormean.cur, col="red",pch=20)

boxplot(importance.individual[1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 1 - 5,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[2,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 5,001 - 10,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[3,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 10,001 - 15,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
boxplot(importance.individual[repeats+1,1:(frac.par.update*nb.k),indices], ylim=c(0,1), names=indices,
        main="Parallel updating \n Iterations 100,001 - 290,000", ylab="Emp. inclusion freq.", xlab="Index of covariate", cex.axis = size_axis, xaxt="n")
axis(1, at = c(1, seq(10, 100, by = 10)), cex.axis = size_axis)
axis(1, at = 1:100, labels=rep("",100))
#points(results$parameter.list[[1]]$priormean.cur, col="red",pch=20)
###################################################################################################







