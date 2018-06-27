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

#Load butachlor concentration values from NAWQA and sample range for r0 simulation
load("Agrochemical_Review/Sims/Range/Butachlor/but_range.RData")

#Load butachlor r0 simulations
load("Agrochemical_Review/Sims/Range/Butachlor/but_r0_sims.RData")

#GGplot theme for manuscripts
source("Agrochemical_Review/Sims/ggplot_theme.R")

#Boxplot of butachlor values from NAWQA
but_box <- data.frame(Chem = "Butachlor",
           Concentration = log(but_vals)) %>% 
  ggplot(aes(x = Chem, y = Concentration)) + 
  geom_boxplot(outlier.shape = 1) + geom_hline(yintercept = log(but_eec), lty = 2) +
  theme_ms() + 
  theme(panel.border = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  coord_flip() + scale_y_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 20)),
                                    labels = c(0.001, 0.01, 0.1, 1, 20),
                                    limits = log(c(0.0001, 20))) +
  ylab("Butachlor Concentration")

but_box

#R0 over butachlor concentration plot
but_r0 <- but_sims %>% 
  ggplot(aes(x = log(but), y = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)))) + 
  geom_line() + geom_vline(xintercept = log(but_eec), lty = 2) +
  theme_ms() + ylab(expression(R[0])) +
  geom_ribbon(aes(x = log(but), 
                  ymin = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)), 
                  ymax = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))), alpha = 0.4) + 
  scale_x_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 20)),
                     labels = c(0.001, 0.01, 0.1, 1, 20),
                     limits = log(c(0.0001, 20))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

but_r0

#Put plots together (code adapted from http://felixfan.github.io/stacking-plots-same-x/)
but1 <- ggplotGrob(but_r0)
but2 <- ggplotGrob(but_box)
# stack the two plots
but <- rbind(but1, but2, size="first") 
# use the largest widths
but$widths <- unit.pmax(but1$widths, but2$widths) 
#Adjust the heights
panels <- but$layout$t[grep("panel", but$layout$name)]
but$heights[panels[1]] <- unit(5, "null") 
but$heights[panels[2]] <- unit(1,"null") 

#Print the plot
grid.draw(but)  
