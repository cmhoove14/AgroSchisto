require(rootSolve)
require(stringr)
require(tidyverse)

source("Agrochemical_Review/Models/r0_functions.R")
source("Agrochemical_Review/Models/model_helper_functions.R")
source("Agrochemical_Review/Models/initial_parameters.R")

test_Ws <- c(seq(0.0001, 1, 0.0001), seq(1, 100, 1))

test_ks <- c(0.001, 0.01, 0.1, 1, 10)

test_phis <- as.data.frame(expand.grid("W" = test_Ws, 
                                       "k" = test_ks)) %>% 
  mutate(phi = map2_dbl(W, k, phi_Wk))

test_phis %>% 
  mutate(k = as.factor(k)) %>% 
  ggplot(aes(x = W, y = phi, col = k)) +
    geom_line(size = 1.2) +
    theme_classic() +
    theme(legend.position = c(0.5,0.5),
          axis.title = element_text(size = 14)) +
    labs(x = "Mean worm burden (W)",
         y = expression(paste("Mating probability (", Phi,")")))
    
    
