library(gridExtra)
library(zoo)
library(tidyverse)

#GGplot theme for manuscripts
source("Agrochemical_Review/Sims/ggplot_theme.R")

#load R0 function
source("Agrochemical_Review/Models/r0_of_q.R")

#Load profenofos r0 simulations and nawqa data
load("Agrochemical_Review/Sims/Range/Profenofos/prof_range.RData")
load("Agrochemical_Review/Sims/Range/Profenofos/prof_r0_sims.RData")
  prof_sims <- prof_sims %>% mutate(Insecticide = "Profenofos",
                                    EEC = prof / prof_eec,
                                    r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                    r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                    r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = prof) %>% 
    filter(EEC <= 1)

#Load malathion r0 simulations
load("Agrochemical_Review/Sims/Range/Malathion/mal_range.RData")
load("Agrochemical_Review/Sims/Range/Malathion/mal_r0_sims.RData")
  mal_sims <- mal_sims %>% mutate(Insecticide = "Malathion",
                                  EEC = mal / mal_eec,
                                  r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                  r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = mal) %>% 
    filter(EEC <= 1.0001)

#Load chlorpyrifos r0 simulations
load("Agrochemical_Review/Sims/Range/Chlorpyrifos/chlor_range.RData")
load("Agrochemical_Review/Sims/Range/Chlorpyrifos/chlor_r0_sims.RData")
  chlor_sims <- chlor_sims %>% mutate(Insecticide = "Chlorpyrifos",
                                      EEC = chlor / chlor_eec,
                                      r0_smooth = rollmean(r0_med, k = 10, align = "left", fill = c(NA,NA,NA)),
                                      r025_smooth = rollmean(r0_25, k = 10, align = "left", fill = c(NA,NA,NA)),
                                      r075_smooth = rollmean(r0_75, k = 10, align = "left", fill = c(NA,NA,NA))) %>% 
    rename(conc = chlor) %>% 
    filter(EEC <= 1)

# All insecticides together
ins_plot <- rbind(prof_sims, mal_sims, chlor_sims) %>% 
  ggplot(aes(x = EEC, y = r0_smooth, col = Insecticide)) + 
  geom_line(size = 1) + #geom_vline(xintercept = log(prof_eec), lty = 2) +
  geom_ribbon(aes(x = EEC, ymin = r025_smooth, ymax = r075_smooth, fill = Insecticide), alpha = 0.25) + 
  theme_ms() + ylab(expression(R[0])) + xlab("Agrochemical Concentration Normalized to Peak EEC") +
  ggtitle(expression(paste("Combined effects of key insecticides on estimates of ", R[0]))) +
  scale_x_continuous(breaks = seq(0,1,0.5),
                     labels = c("0", "0.5 Peak EEC", "Peak EEC"),
                     limits = c(0, 1.0001)) +
  scale_y_continuous(breaks = seq(1.25, 3.75, .5),
                     labels = seq(1.25, 3.75, .5),
                     limits = c(1, 3.75)) +
  geom_hline(yintercept = r0.fix()[3], lty = 2)

ins_plot