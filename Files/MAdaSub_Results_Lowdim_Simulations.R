load("./Files/R_simulation_saves/results_lowdim_03.RData")
load("./Files/R_simulation_saves/results_lowdim_09.RData")
load("./Files/R_simulation_saves/results_lowdim_099.RData")

Iterations <- 20000

n_sim <- length(results_lowdim_03$acceptance_rates_MAdaSub)

#methodnames = c("MA (a)", "MA (b)", "MC3", "MC3 swap", "wTGS")
#methodnames = c("MA (a)", "MA (b)", "MC3", "MC3 swap", "TGS", "RB-TGS", "wTGS", "RB-wTGS")
methodnames = c("MA (a)", "MA (b)", "MC3", "MC3 swap", "wTGS", "wTGS R-B")


font = 1
par(las=1,lwd=1,cex.main=1.2, cex.lab=1, cex.axis=1.3)
options(scipen=5)

#dev.off()
#win.graph(width = 14, height = 8, pointsize = 8)

########################## Revised plot

library(ggplot2)
library(RColorBrewer)

values <- c(results_lowdim_03$acceptance_rates_MAdaSub, 
                      results_lowdim_03$acceptance_rates_PO, 
                      results_lowdim_03$acceptance_rates_MC3, 
                      results_lowdim_03$acceptance_rates_MC3_swap, 
                      #results_lowdim_03$acceptance_rates_TGS,
                      #results_lowdim_03$acceptance_rates_TGS,
                      results_lowdim_03$acceptance_rates_wTGS,
                      results_lowdim_03$acceptance_rates_wTGS,
                      results_lowdim_09$acceptance_rates_MAdaSub, 
                      results_lowdim_09$acceptance_rates_PO, 
                      results_lowdim_09$acceptance_rates_MC3, 
                      results_lowdim_09$acceptance_rates_MC3_swap, 
                      #results_lowdim_09$acceptance_rates_TGS,
                      #results_lowdim_09$acceptance_rates_TGS,
                      results_lowdim_09$acceptance_rates_wTGS,
                      results_lowdim_09$acceptance_rates_wTGS,
                      results_lowdim_099$acceptance_rates_MAdaSub, 
                      results_lowdim_099$acceptance_rates_PO, 
                      results_lowdim_099$acceptance_rates_MC3, 
                      results_lowdim_099$acceptance_rates_MC3_swap, 
                      #results_lowdim_099$acceptance_rates_TGS,
                      #results_lowdim_099$acceptance_rates_TGS,
                      results_lowdim_099$acceptance_rates_wTGS,
                      results_lowdim_099$acceptance_rates_wTGS,
                          results_lowdim_03$Iterations_converged_MAdaSub, 
                          results_lowdim_03$Iterations_converged_PO, 
                          results_lowdim_03$Iterations_converged_MC3, 
                          results_lowdim_03$Iterations_converged_MC3_swap, 
                       #   results_lowdim_03$Iterations_converged_TGS_vanilla,
                       #   results_lowdim_03$Iterations_converged_TGS,
                          results_lowdim_03$Iterations_converged_wTGS_vanilla,
                          results_lowdim_03$Iterations_converged_wTGS,
                          results_lowdim_09$Iterations_converged_MAdaSub, 
                          results_lowdim_09$Iterations_converged_PO, 
                          results_lowdim_09$Iterations_converged_MC3, 
                          results_lowdim_09$Iterations_converged_MC3_swap, 
                      #    results_lowdim_09$Iterations_converged_TGS_vanilla,
                      #    results_lowdim_09$Iterations_converged_TGS,
                          results_lowdim_09$Iterations_converged_wTGS_vanilla,
                          results_lowdim_09$Iterations_converged_wTGS,
                          results_lowdim_099$Iterations_converged_MAdaSub, 
                          results_lowdim_099$Iterations_converged_PO, 
                          results_lowdim_099$Iterations_converged_MC3, 
                          results_lowdim_099$Iterations_converged_MC3_swap, 
                      #    results_lowdim_099$Iterations_converged_TGS_vanilla,
                      #    results_lowdim_099$Iterations_converged_TGS,
                          results_lowdim_099$Iterations_converged_wTGS_vanilla,
                          results_lowdim_099$Iterations_converged_wTGS)

methods <- rep(methodnames, times = 3*2, each = n_sim)
Correlation <- rep(c("0.3", "0.9", "0.99"), times = 2, each = n_sim*length(methodnames))
output <- rep(c("Acceptance rate", "Iterations for PIP convergence"), each = n_sim*length(methodnames)*3)

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







