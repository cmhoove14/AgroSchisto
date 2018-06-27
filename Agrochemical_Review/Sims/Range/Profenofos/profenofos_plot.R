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

#Load profenofos concentration values from NAWQA and sample range for r0 simulation
load("Agrochemical_Review/Sims/Range/Profenofos/prof_range.RData")

#Load profenofos r0 simulations
load("Agrochemical_Review/Sims/Range/Profenofos/prof_r0_sims.RData")

#GGplot theme for manuscripts
source("Agrochemical_Review/Sims/ggplot_theme.R")

#Boxplot of profenofos values from NAWQA
prof_box <- data.frame(Chem = "Profenofos",
           Concentration = log(prof_vals)) %>% 
  ggplot(aes(x = Chem, y = Concentration)) + geom_boxplot(outlier.shape = 1) + 
  theme_ms() + geom_hline(yintercept = log(prof_eec), lty = 2) +
  theme(panel.border = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  coord_flip() + scale_y_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 10, 20)),
                                    labels = c(0.001, 0.01, 0.1, 1, 10, 20),
                                    limits = log(c(0.001, 20))) +
  ylab(expression(paste("Profenofos Concentration (", mu, "g/L)")))

prof_box

#R0 over profenofos concentration plot
prof_r0 <- prof_sims %>% 
  ggplot(aes(x = log(prof), y = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)))) + 
  geom_line() + geom_vline(xintercept = log(prof_eec), lty = 2) +
  theme_ms() + ylab(expression(R[0])) +
  geom_ribbon(aes(x = log(prof), 
                  ymin = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)), 
                  ymax = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))), alpha = 0.4) + 
  scale_x_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 10, 20)),
                     labels = c(0.001, 0.01, 0.1, 1, 10, 20),
                     limits = log(c(0.001, 20))) +
  scale_y_continuous(breaks = seq(1.5, 3.75, 0.5),
                     labels = seq(1.5, 3.75, 0.5),
                     limits = c(1.5, 3.75)) +
  geom_hline(yintercept = prof_sims$r0_med[1], lty = 3) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

prof_r0

#Put plots together (code adapted from http://felixfan.github.io/stacking-plots-same-x/)
prof1 <- ggplotGrob(prof_r0)
prof2 <- ggplotGrob(prof_box)
# stack the two plots
prof <- rbind(prof1, prof2, size="first") 
# use the largest widths
prof$widths <- unit.pmax(prof1$widths, prof2$widths) 
#Adjust the heights
panels <- prof$layout$t[grep("panel", prof$layout$name)]
prof$heights[panels[1]] <- unit(5, "null") 
prof$heights[panels[2]] <- unit(1,"null") 

#Print the plot
grid.draw(prof)  
