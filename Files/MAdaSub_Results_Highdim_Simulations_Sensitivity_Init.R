
nb.k = 5 # 5 
n = 500 # n = 500 or n = 1000
p = 500 # p = 500 or p = 5000
SNR = 0.5
repeats = 50 # 50
if (p==500) {
  M = 1000
}
if (p==5000) {
  M = 10000
}

epsilon = 1/p
#load(paste("MAdaSub_results_highdim_sim_n", n, "_p", p, "_SNR", SNR, "_CPU", nb.k, "_Rep", repeats, "_M", M, ".RData", sep=""))
load(paste("./Files/R_simulation_saves/MAdaSub_results_highdim_parallel_without2_sim_n", n, "_p", p, "_SNR", SNR, "_CPU", nb.k, "_Rep", repeats, "_M", M, ".RData", sep=""))


var_MAdaSub_parallel <- apply(results$Estimates_MAdaSub_parallel, 2, var)
var_MAdaSub_parallel_without <- apply(results$Estimates_MAdaSub_parallel_without, 2, var)
var_MAdaSub_parallel_only_rj <- apply(results$Estimates_MAdaSub_parallel_only_rj, 2, var)
var_MC3 <- apply(results$Estimates_MC3, 2, var)

r_hat_parallel <- var_MC3 * median(results$Times_MC3) / (var_MAdaSub_parallel * median(results$Times_MAdaSub_parallel))
r_hat_parallel_without <- var_MC3 * median(results$Times_MC3) / (var_MAdaSub_parallel_without * median(results$Times_MAdaSub_parallel_without))
r_hat_parallel_only_rj <- var_MC3 * median(results$Times_MC3) / (var_MAdaSub_parallel_only_rj * median(results$Times_MAdaSub_parallel_only_rj))

n_best <- 20
highest_PIP <- order(colMeans(rbind(results$Estimates_MAdaSub_parallel,
                                    results$Estimates_MAdaSub_parallel_without,
                                    results$Estimates_MAdaSub_parallel_only_rj,
                                    results$Estimates_MC3)), decreasing=TRUE)[1:n_best]

sorted_PIP <- sort(colMeans(rbind(results$Estimates_MAdaSub_parallel,
                                  results$Estimates_MAdaSub_parallel_without,
                                  results$Estimates_MAdaSub_parallel_only_rj,
                                  results$Estimates_MC3)), decreasing=TRUE)
#sorted_PIP[21]*100
#median(sorted_PIP)*100

median(r_hat_parallel_without[highest_PIP], na.rm=TRUE)
median(r_hat_parallel_only_rj[highest_PIP], na.rm=TRUE)
median(r_hat_parallel[highest_PIP], na.rm=TRUE)

median(results$Accrate_MAdaSub_parallel_without)
median(results$Accrate_MAdaSub_parallel_only_rj)
median(results$Accrate_MAdaSub_parallel)

median(results$Accrate_MC3)

median(r_hat_parallel_without, na.rm=TRUE)
median(r_hat_parallel_only_rj, na.rm=TRUE)
median(r_hat_parallel, na.rm=TRUE)





