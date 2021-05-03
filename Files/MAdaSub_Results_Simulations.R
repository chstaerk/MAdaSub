
load("./Files/R_simulation_saves/results_lowdim_03_epsilon.RData")
load("./Files/R_simulation_saves/results_lowdim_09_epsilon.RData")
load("./Files/R_simulation_saves/results_lowdim_099_epsilon.RData")

methodnames = c("MAdaSub (a)", "MAdaSub (b)", "MC3")


font = 1
par(las=1,lwd=1,cex.main=font, cex.lab=font, cex.axis=1.3)
par(cex.main=1.2)
options(scipen=5)

par(mfcol=c(2,3))
par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))
results = results_lowdim_03
#expression(paste("Acceptance rates", rho," = 0.3"))
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, main = expression(bold(atop("Acceptance rate", rho ~" = 0.3"))), ylim = c(0,1), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.3"))), ylim = c(0,Iterations+1), names = methodnames)

results = results_lowdim_09
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, main = expression(bold(atop("Acceptance rate", rho ~" = 0.9"))), ylim = c(0,1), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.9"))), ylim = c(0,Iterations+1), names = methodnames)

results = results_lowdim_099
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, main = expression(bold(atop("Acceptance rate", rho ~" = 0.99"))), ylim = c(0,1), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.99"))), ylim = c(0,Iterations+1), names = methodnames)
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))




#####################################################################

font = 1
par(las=1,lwd=1,cex.main=font, cex.lab=font, cex.axis=1.3)
par(cex.main=1.2)
options(scipen=5)

par(mfcol=c(3,3))
par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))
results = results_lowdim_03
#expression(paste("Acceptance rates", rho," = 0.3"))
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, main = expression(bold(atop("Acceptance rate", rho ~" = 0.3"))), ylim = c(0,1), names = methodnames)
boxplot(results$effective_sizes_median_MAdaSub, results$effective_sizes_median_PO, results$effective_sizes_median_MC3, main = expression(bold(atop("Median efective sample size", rho ~" = 0.3"))), ylim = c(0,Iterations), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.3"))), ylim = c(0,Iterations+1), names = methodnames)

results = results_lowdim_09
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, main = expression(bold(atop("Acceptance rate", rho ~" = 0.9"))), ylim = c(0,1), names = methodnames)
boxplot(results$effective_sizes_median_MAdaSub, results$effective_sizes_median_PO, results$effective_sizes_median_MC3, main = expression(bold(atop("Median efective sample size", rho ~" = 0.9"))), ylim = c(0,Iterations), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.9"))), ylim = c(0,Iterations+1), names = methodnames)

results = results_lowdim_099
boxplot(results$acceptance_rates_MAdaSub, results$acceptance_rates_PO, results$acceptance_rates_MC3, main = expression(bold(atop("Acceptance rate", rho ~" = 0.99"))), ylim = c(0,1), names = methodnames)
boxplot(results$effective_sizes_median_MAdaSub, results$effective_sizes_median_PO, results$effective_sizes_median_MC3, main = expression(bold(atop("Median efective sample size", rho ~" = 0.99"))), ylim = c(0,Iterations), names = methodnames)
boxplot(results$Iterations_converged_MAdaSub, results$Iterations_converged_PO, results$Iterations_converged_MC3, main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.99"))), ylim = c(0,Iterations+1), names = methodnames)
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

