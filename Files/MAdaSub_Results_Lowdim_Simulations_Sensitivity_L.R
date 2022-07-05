library("RColorBrewer")

load("./Files/R_simulation_saves/results_lowdim_03_Lj.RData")
load("./Files/R_simulation_saves/results_lowdim_09_Lj.RData")
load("./Files/R_simulation_saves/results_lowdim_099_Lj.RData")

Iterations = 20000
p = 20
Ljs = c(p/20, p/10, p/5, p/2, p, 2*p, 5*p, 10*p, 20*p)

## For L=p use the same result from original run (differences only due to slightly varying computational time/ random seeds)
load("results_lowdim_03.RData")
boxplot(results_lowdim_03_Lj$acceptance_rates_MAdaSub[,5],  results_lowdim_03$acceptance_rates_MAdaSub)
wilcox.test(results_lowdim_03_Lj$acceptance_rates_MAdaSub[,5], results_lowdim_03$acceptance_rates_MAdaSub)
boxplot(results_lowdim_03_Lj$Iterations_converged_MAdaSub[,5],  results_lowdim_03$Iterations_converged_MAdaSub)
wilcox.test(results_lowdim_03_Lj$Iterations_converged_MAdaSub[,5],  results_lowdim_03$Iterations_converged_MAdaSub)
results_lowdim_03_Lj$acceptance_rates_MAdaSub[,5] <- results_lowdim_03$acceptance_rates_MAdaSub
results_lowdim_03_Lj$Iterations_converged_MAdaSub[,5] <- results_lowdim_03$Iterations_converged_MAdaSub

load("results_lowdim_09.RData")
boxplot(results_lowdim_09_Lj$acceptance_rates_MAdaSub[,5],  results_lowdim_09$acceptance_rates_MAdaSub)
wilcox.test(results_lowdim_09_Lj$acceptance_rates_MAdaSub[,5],  results_lowdim_09$acceptance_rates_MAdaSub)
boxplot(results_lowdim_09_Lj$Iterations_converged_MAdaSub[,5],  results_lowdim_09$Iterations_converged_MAdaSub)
wilcox.test(results_lowdim_09_Lj$Iterations_converged_MAdaSub[,5],  results_lowdim_09$Iterations_converged_MAdaSub)
results_lowdim_09_Lj$acceptance_rates_MAdaSub[,5] <- results_lowdim_09$acceptance_rates_MAdaSub
results_lowdim_09_Lj$Iterations_converged_MAdaSub[,5] <- results_lowdim_09$Iterations_converged_MAdaSub

load("results_lowdim_099.RData")
boxplot(results_lowdim_099_Lj$acceptance_rates_MAdaSub[,5],  results_lowdim_099$acceptance_rates_MAdaSub)
wilcox.test(results_lowdim_099_Lj$acceptance_rates_MAdaSub[,5],  results_lowdim_099$acceptance_rates_MAdaSub)
boxplot(results_lowdim_099_Lj$Iterations_converged_MAdaSub[,5],  results_lowdim_099$Iterations_converged_MAdaSub)
wilcox.test(results_lowdim_099_Lj$Iterations_converged_MAdaSub[,5],  results_lowdim_099$Iterations_converged_MAdaSub)
results_lowdim_099_Lj$acceptance_rates_MAdaSub[,5] <- results_lowdim_099$acceptance_rates_MAdaSub
results_lowdim_099_Lj$Iterations_converged_MAdaSub[,5] <- results_lowdim_099$Iterations_converged_MAdaSub
####################################################################

font=1.3
par(las=1,lwd=1,cex.main=font, cex.lab=font, cex.axis=1.18)
par(mar=c(5, 4, 4, 1), mgp = c(2.5, 1, 0))

par(mfcol=c(2,3))
boxplot(results_lowdim_03_Lj$acceptance_rates_MAdaSub, names=Ljs, ylim=c(0,1), main = expression(bold(atop("Acceptance rate", rho ~" = 0.3"))), xlab="L",
        col=brewer.pal(n = length(Ljs), name = "RdYlBu"))
boxplot(results_lowdim_03_Lj$Iterations_converged_MAdaSub, names=Ljs, ylim = c(0,Iterations+1), main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.3"))), xlab="L",
        col=brewer.pal(n = length(Ljs), name = "RdYlBu"))

boxplot(results_lowdim_09_Lj$acceptance_rates_MAdaSub, names=Ljs, ylim=c(0,1), main = expression(bold(atop("Acceptance rate", rho ~" = 0.9"))), xlab="L",
        col=brewer.pal(n = length(Ljs), name = "RdYlBu"))
boxplot(results_lowdim_09_Lj$Iterations_converged_MAdaSub, names=Ljs, ylim = c(0,Iterations+1), main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.9"))), xlab="L",
        col=brewer.pal(n = length(Ljs), name = "RdYlBu"))

boxplot(results_lowdim_099_Lj$acceptance_rates_MAdaSub, names=Ljs, ylim=c(0,1), main = expression(bold(atop("Acceptance rate", rho ~" = 0.99"))), xlab="L",
        col=brewer.pal(n = length(Ljs), name = "RdYlBu"))
boxplot(results_lowdim_099_Lj$Iterations_converged_MAdaSub, names=Ljs, ylim = c(0,Iterations+1), main = expression(bold(atop("Iterations for PIP convergence", rho ~" = 0.99"))), xlab="L",
        col=brewer.pal(n = length(Ljs), name = "RdYlBu"))

par(mar=c(4, 4, 4, 1), mgp = c(2.5, 1, 0))

library(ggplot2)
library(RColorBrewer)

values <- c(c(results_lowdim_03_Lj$acceptance_rates_MAdaSub), 
            c(results_lowdim_09_Lj$acceptance_rates_MAdaSub), 
            c(results_lowdim_099_Lj$acceptance_rates_MAdaSub),
            c(results_lowdim_03_Lj$Iterations_converged_MAdaSub), 
            c(results_lowdim_09_Lj$Iterations_converged_MAdaSub), 
            c(results_lowdim_099_Lj$Iterations_converged_MAdaSub))

n_sim <- length(results_lowdim_03_Lj$acceptance_rates_MAdaSub[,1])

methods <- as.factor(paste("L =", rep(Ljs, times = 3*2, each = n_sim)))
Correlation <- rep(c("0.3", "0.9", "0.99"), times = 2, each = n_sim*length(Ljs))
output <- rep(c("Acceptance rate", "Iterations for PIP convergence"), each = n_sim*length(Ljs)*3)

dat_results <- data.frame(methods, Correlation, output, values)

dim(dat_results)

ggplot(dat_results,                          
       aes(x = methods,
           y = values,
           fill = Correlation)) +
  geom_boxplot() +
  facet_wrap(~output, ncol=1, scales = "free") + 
  labs(fill = expression( "Correlation" ~ rho )) +
  expand_limits(x = 0, y = 0) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_fill_brewer(palette="YlGnBu") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x=element_blank()) +
  theme(strip.background =element_rect(fill="white")) +
  theme(panel.spacing = unit(1.5, "lines")) +
  theme( axis.text = element_text( size = 15 ),
         axis.text.x = element_text( size = 15 ),
         strip.text = element_text(size = 18),
         legend.text =  element_text(size = 15),
         legend.title = element_text(size=15))



