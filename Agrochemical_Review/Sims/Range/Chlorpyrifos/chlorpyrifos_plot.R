#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

library(gridExtra)
library(zoo)
library(tidyverse)

#Load chlorpyrifos concentration values from NAWQA and sample range for r0 simulation
load("Agrochemical_Review/Sims/Range/Chlorpyrifos/chlor_range.RData")

#Load chlorpyrifos r0 simulations
load("Agrochemical_Review/Sims/Range/Chlorpyrifos/chlor_r0_sims.RData")

#GGplot theme for manuscripts
source("Agrochemical_Review/Sims/ggplot_theme.R")

#Boxplot of chlorpyrifos values from NAWQA
chlor_box <- data.frame(Chem = "Chlorpyrifos",
           Concentration = log(chlor_vals)) %>% 
  ggplot(aes(x = Chem, y = Concentration)) + geom_boxplot(outlier.shape = 1) + 
  theme_ms() + geom_hline(yintercept = log(chlor_eec), lty = 2) +
  theme(panel.border = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  coord_flip() + scale_y_continuous(breaks = log(c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100)),
                                    labels = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100),
                                    limits = log(c(0.0001, 100))) +
  ylab(expression(paste("Chlorpyrifos Concentration (", mu, "g/L)")))

chlor_box

#R0 over chlorpyrifos concentration plot
chlor_r0 <- chlor_sims %>% 
  ggplot(aes(x = log(chlor), y = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)))) + 
  geom_line() + geom_vline(xintercept = log(chlor_eec), lty = 2) +
  theme_ms() + ylab(expression(R[0])) +
  geom_ribbon(aes(x = log(chlor), 
                  ymin = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)), 
                  ymax = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))), alpha = 0.4) + 
  scale_x_continuous(breaks = log(c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100)),
                     labels = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100),
                     limits = log(c(0.0001, 100))) +
  scale_y_continuous(breaks = c(1.25, 1.75, 2.25, 2.75, 3.25),
                     labels = c(1.25, 1.75, 2.25, 2.75, 3.25),
                     limits = (c(1.25, 3.25))) +
  geom_hline(yintercept = chlor_sims$r0_med[1], lty = 3) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

chlor_r0

#Put plots together (code adapted from http://felixfan.github.io/stacking-plots-same-x/)
chlor1 <- ggplotGrob(chlor_r0)
chlor2 <- ggplotGrob(chlor_box)
# stack the two plots
chlor <- rbind(chlor1, chlor2, size="first") 
# use the largest widths
chlor$widths <- unit.pmax(chlor1$widths, chlor2$widths) 
#Adjust the heights
panels <- chlor$layout$t[grep("panel", chlor$layout$name)]
chlor$heights[panels[1]] <- unit(5, "null") 
chlor$heights[panels[2]] <- unit(1,"null") 

#Print the plot
grid.draw(chlor)  
