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
  atr_sims <- atr_sims %>% mutate(Agrochemical = "Atrazine",
                                  Class = "Herbicide",
                                  EEC = atr / atr_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = atr) %>% 
    filter(EEC <= 1.2)

#Load butachlor r0 simulations and nawqa data
load("Agrochemical_Review/Sims/Range/Butachlor/but_r0_sims.RData")
load("Agrochemical_Review/Sims/Range/Butachlor/but_range.RData")
  but_sims <- but_sims %>% mutate(Agrochemical = "Butachlor",
                                  Class = "Herbicide",
                                  EEC = but / but_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = but) %>% 
    filter(EEC <= 1.2)

#Load butralin r0 simulations and nawqa data
load("Agrochemical_Review/Sims/Range/Butralin/btr_r0_sims.RData")
load("Agrochemical_Review/Sims/Range/Butralin/btr_range.RData")
  btr_sims <- btr_sims %>% mutate(Agrochemical = "Butralin",
                                  Class = "Herbicide",
                                  EEC = btr / btr_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = btr) %>% 
    filter(EEC <= 1.2)
 
#Load glyphosate r0 simulations and nawqa data
load("Agrochemical_Review/Sims/Range/Glyphosate/gly_r0_sims.RData")
load("Agrochemical_Review/Sims/Range/Glyphosate/gly_range.RData")
  gly_sims <- gly_sims %>% mutate(Agrochemical = "Glyphosate",
                                  Class = "Herbicide",
                                  EEC = gly / gly_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = gly) %>% 
    filter(EEC <= 1.2)

#Load profenofos r0 simulations and nawqa data
load("Agrochemical_Review/Sims/Range/Profenofos/prof_range.RData")
load("Agrochemical_Review/Sims/Range/Profenofos/prof_r0_sims.RData")
  prof_sims <- prof_sims %>% mutate(Agrochemical = "Profenofos",
                                    Class = "Insecticide",
                                    EEC = prof / prof_eec,
                                    r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                    r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                    r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = prof) %>% 
    filter(EEC <= 1)

#Load malathion r0 simulations
load("Agrochemical_Review/Sims/Range/Malathion/mal_range.RData")
load("Agrochemical_Review/Sims/Range/Malathion/mal_r0_sims.RData")
  mal_sims <- mal_sims %>% mutate(Agrochemical = "Malathion",
                                  Class = "Insecticide",
                                  EEC = mal / mal_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = mal) %>% 
    filter(EEC <= 1.0001)

#Load chlorpyrifos r0 simulations
load("Agrochemical_Review/Sims/Range/Chlorpyrifos/chlor_range.RData")
load("Agrochemical_Review/Sims/Range/Chlorpyrifos/chlor_r0_sims.RData")
  chlor_sims <- chlor_sims %>% mutate(Agrochemical = "Chlorpyrifos",
                                      Class = "Insecticide",
                                      EEC = chlor / chlor_eec,
                                      r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                      r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                      r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = chlor) %>% 
    filter(EEC <= 1)
     
# All herbicides together
agroc_plot <- rbind(chlor_sims, prof_sims, mal_sims, 
                    atr_sims, but_sims, btr_sims, gly_sims) %>% 
  ggplot(aes(x = EEC, y = r0_smooth, col = Agrochemical)) + 
  geom_line(size = 1) + facet_grid(Class ~ .) +
  geom_ribbon(aes(x = EEC, ymin = r025_smooth, ymax = r075_smooth, fill = Agrochemical), alpha = 0.25) + 
  theme_ms() +
  ylab(expression(R[0])) + xlab("Agrochemical Concentration Normalized to Peak EEC") +
  ggtitle(expression(paste("Combined effects of key agrochemicals on estimates of ", R[0]))) +
  scale_x_continuous(breaks = seq(0,1,0.5),
                     labels = c("0", "0.5 Peak EEC", "Peak EEC"),
                     limits = c(0, 1.001)) +
  geom_hline(yintercept = r0.fix()[3], lty = 2)
  
agroc_plot
