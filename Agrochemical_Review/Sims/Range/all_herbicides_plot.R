library(gridExtra)
library(zoo)
library(tidyverse)

#GGplot theme for manuscripts
source("Agrochemical_Review/Sims/ggplot_theme.R")

#load R0 function
source("Agrochemical_Review/Models/r0_of_q.R")

#Load atrazine r0 simulations and nawqa data
load("Agrochemical_Review/Sims/Range/Atrazine/atr_r0_sims.RData")
load("Agrochemical_Review/Sims/Range/Atrazine/atr_range.RData")
  atr_sims <- atr_sims %>% mutate(Herbicide = "Atrazine",
                                  EEC = atr / atr_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = atr) %>% 
    filter(EEC <= 1.2)

#Load butachlor r0 simulations and nawqa data
load("Agrochemical_Review/Sims/Range/Butachlor/but_r0_sims.RData")
load("Agrochemical_Review/Sims/Range/Butachlor/but_range.RData")
  but_sims <- but_sims %>% mutate(Herbicide = "Butachlor",
                                  EEC = but / but_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = but) %>% 
    filter(EEC <= 1.2)

#Load butralin r0 simulations and nawqa data
load("Agrochemical_Review/Sims/Range/Butralin/btr_r0_sims.RData")
load("Agrochemical_Review/Sims/Range/Butralin/btr_range.RData")
  btr_sims <- btr_sims %>% mutate(Herbicide = "Butralin",
                                  EEC = btr / btr_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = btr) %>% 
    filter(EEC <= 1.2)
 
#Load glyphosate r0 simulations and nawqa data
load("Agrochemical_Review/Sims/Range/Glyphosate/gly_r0_sims.RData")
load("Agrochemical_Review/Sims/Range/Glyphosate/gly_range.RData")
  gly_sims <- gly_sims %>% mutate(Herbicide = "Glyphosate",
                                  EEC = gly / gly_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = gly) %>% 
    filter(EEC <= 1.2)
   
# All herbicides together
herb_plot <- rbind(atr_sims, but_sims, btr_sims, gly_sims) %>% 
  ggplot(aes(x = EEC, y = r0_smooth, col = Herbicide)) + 
  geom_line(size = 1) + #geom_vline(xintercept = log(prof_eec), lty = 2) +
  geom_ribbon(aes(x = EEC, ymin = r025_smooth, ymax = r075_smooth, fill = Herbicide), alpha = 0.25) + 
  theme_ms() + ylab(expression(R[0])) + xlab("Agrochemical Concentration Normalized to Peak EEC") +
  ggtitle(expression(paste("Combined effects of key herbicides on estimates of ", R[0]))) +
  scale_x_continuous(breaks = seq(0,1,0.5),
                     labels = c("0", "0.5 Peak EEC", "Peak EEC"),
                     limits = c(0, 1.001)) +
  scale_y_continuous(breaks = seq(0, 2.5, .5),
                     labels = seq(0, 2.5, .5),
                     limits = c(-0.001, 2.5)) +
  geom_hline(yintercept = r0.fix()[3], lty = 2)
  
herb_plot