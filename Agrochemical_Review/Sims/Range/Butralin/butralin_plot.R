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

#Load butralin concentration values from NAWQA and sample range for r0 simulation
load("Agrochemical_Review/Sims/Range/Butralin/btr_range.RData")

#Load butralin r0 simulations
load("Agrochemical_Review/Sims/Range/Butralin/btr_r0_sims.RData")

#GGplot theme for manuscripts
source("Agrochemical_Review/Sims/ggplot_theme.R")

#Boxplot of butralin values from NAWQA
btr_box <- data.frame(Chem = "Butralin",
           Concentration = log(btr_vals)) %>% 
  ggplot(aes(x = Chem, y = Concentration)) + 
  geom_boxplot(outlier.shape = 1) + geom_hline(yintercept = log(btr_eec), lty = 2) +
  theme_ms() + 
  theme(panel.border = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  coord_flip() + scale_y_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 20)),
                                    labels = c(0.001, 0.01, 0.1, 1, 20),
                                    limits = log(c(0.0001, 20))) +
  ylab(expression(paste("Butralin Concentration (", mu, "g/L)")))

btr_box

#R0 over butralin concentration plot
btr_r0 <- btr_sims %>% 
  ggplot(aes(x = log(btr), y = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)))) + 
  geom_line() + geom_vline(xintercept = log(btr_eec), lty = 2) +
  theme_ms() + ylab(expression(R[0])) +
  geom_ribbon(aes(x = log(btr), 
                  ymin = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)), 
                  ymax = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))), alpha = 0.4) + 
  scale_x_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 20)),
                     labels = c(0.001, 0.01, 0.1, 1, 20),
                     limits = log(c(0.0001, 20))) +
  scale_y_continuous(breaks = c(1.5, 1.75, 2.0, 2.25),
                     labels = c(1.5, 1.75, 2.0, 2.25),
                     limits = (c(1.5, 2.25))) +
  geom_hline(yintercept = btr_sims$r0_med[1], lty = 3) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

btr_r0

#Put plots together (code adapted from http://felixfan.github.io/stacking-plots-same-x/)
btr1 <- ggplotGrob(btr_r0)
btr2 <- ggplotGrob(btr_box)
# stack the two plots
btr <- rbind(btr1, btr2, size="first") 
# use the largest widths
btr$widths <- unit.pmax(btr1$widths, btr2$widths) 
#Adjust the heights
panels <- btr$layout$t[grep("panel", btr$layout$name)]
btr$heights[panels[1]] <- unit(5, "null") 
btr$heights[panels[2]] <- unit(1,"null") 

#Print the plot
grid.draw(btr)  
