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

#Load glyphosate concentration values from NAWQA and sample range for r0 simulation
load("Agrochemical_Review/Sims/Range/Glyphosate/gly_range.RData")

#Load glyphosate r0 simulations
load("Agrochemical_Review/Sims/Range/Glyphosate/gly_r0_sims.RData")

#GGplot theme for manuscripts
theme_ms <- function(base_size=12, base_family="Times New Roman") {
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family)+
      theme(text=element_text(color="black"),
            axis.title=element_text(face="bold", size = rel(1.3)),
            axis.text=element_text(size = rel(1), color = "black"),
            legend.title=element_text(face="bold"),
            legend.text=element_text(face="bold"),
            legend.background=element_rect(fill="transparent"),
            legend.key.size = unit(0.8, 'lines'),
            panel.border=element_rect(color="black",size=1)
    ))
}

#Boxplot of glyphosate values from NAWQA
gly_box <- data.frame(Chem = "Glyphosate",
           Concentration = log(gly_vals)) %>% 
  ggplot(aes(x = Chem, y = Concentration)) + geom_boxplot(outlier.shape = 1) + 
  theme_ms() + 
  theme(panel.border = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  coord_flip() + scale_y_continuous(breaks = log(c(0.1, 1, 10, 30)),
                                    labels = c(0.1, 1, 10, 30),
                                    limits = log(c(0.1, 30))) +
  ylab("Glyphosate Concentration")

gly_box

#R0 over glyphosate concentration plot
gly_r0 <- gly_sims %>% 
  ggplot(aes(x = log(gly), y = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)))) + 
  geom_line() + theme_ms() + ylab(expression(R[0])) +
  geom_ribbon(aes(x = log(gly), 
                  ymin = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)), 
                  ymax = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))), alpha = 0.4) + 
  scale_x_continuous(breaks = log(c(0.1, 1, 10, 30)),
                     labels = c(0.1, 1, 10, 30),
                     limits = log(c(0.1, 30))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

gly_r0

#Put plots together (code adapted from http://felixfan.github.io/stacking-plots-same-x/)
gly1 <- ggplotGrob(gly_r0)
gly2 <- ggplotGrob(gly_box)
# stack the two plots
gly <- rbind(gly1, gly2, size="first") 
# use the largest widths
gly$widths <- unit.pmax(gly1$widths, gly2$widths) 
#Adjust the heights
panels <- gly$layout$t[grep("panel", gly$layout$name)]
gly$heights[panels[1]] <- unit(5, "null") 
gly$heights[panels[2]] <- unit(1,"null") 

#Print the plot
grid.draw(gly)  
