
nb.k = 5 # 5 
n = 500
p = 500
SNR = 3
repeats = 50 # 50
M = 1000 # 1000 or 4000 or 10000
epsilon = 1/p
load(paste("./Files/R_simulation_saves/MAdaSub_results_highdim_Ljs_sim_n", n, "_p", p, "_SNR", SNR, "_CPU", nb.k, "_Rep", repeats, "_M", M, ".RData", sep=""))

var_MAdaSub <- apply(results$Estimates_MAdaSub, c(2,3), var)
var_MAdaSub_parallel <- apply(results$Estimates_MAdaSub_parallel, c(2,3), var)
var_MC3 <- apply(results$Estimates_MC3, 2, var)

r_hat <- var_MC3 * median(results$Times_MC3) / (var_MAdaSub %*% diag(apply(results$Times_MAdaSub, 2, median)))
r_hat_parallel <- var_MC3 * median(results$Times_MC3) / (var_MAdaSub_parallel %*% diag(apply(results$Times_MAdaSub_parallel, 2, median)))

#apply(r_hat, 2, median)
#apply(r_hat_parallel, 2, median)

######################################
n_best <- 20
highest_PIP <- order(colMeans(rbind(apply(results$Estimates_MAdaSub_parallel, 2, mean),
                     apply(results$Estimates_MAdaSub, 2, mean),
                     apply(results$Estimates_MC3, 2, mean))), decreasing=TRUE)[1:n_best]
#####################################


apply(r_hat[highest_PIP,], 2, median, na.rm=TRUE)
apply(r_hat_parallel[highest_PIP,], 2, median, na.rm=TRUE)

#apply(r_hat, 2, median, na.rm=TRUE)
#apply(r_hat_parallel, 2, median, na.rm=TRUE)

apply(results$Accrate_MAdaSub, 2, median) * 100
apply(results$Accrate_MAdaSub_parallel, 2, median) * 100
#median(results$Accrate_MC3)
results$Ljs

