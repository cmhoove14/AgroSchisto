require(rootSolve)
require(stringr)
require(tidyverse)

source("Agrochemical_Review/Models/r0_functions.R")
source("Agrochemical_Review/Models/model_helper_functions.R")
source("Agrochemical_Review/Models/initial_parameters.R")

RFxSum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv") %>% 
  arrange(Study) %>% 
  group_by(Study) %>% 
  mutate(ID = row_number(),
         Study_ID_Unique = paste0(Study, "_", ID)) %>% ungroup()

#Create labels field for complex axis labels of forestplot
  max_study_string <- max(str_length(RFxSum$study_long))
  max_chem_string <- max(str_length(RFxSum$Chemical))
  max_species_string <- max(str_length(RFxSum$Species))
  
  RFxSum <- RFxSum %>% 
    mutate(parameters_unicode = case_when(parameter == "mupq" ~ "expression(mu[P])",
                                          parameter == "munq" ~ "expression(mu[N])",
                                          parameter == "Knq" ~ "expression(K[N])",
                                          parameter == "psiq" ~ "expression(psi[P])",
                                          parameter == "picq" ~ "expression(pi[C])",
                                          parameter == "pimq" ~ "expression(pi[M])",
                                          parameter == "fnq" ~ "expression(f[N])",
                                          parameter == "thetaq" ~ "expression(theta)",
                                          parameter == "vq" ~ "expression(v)",
                                          parameter == "fpq" ~ "expression(f[P])") ,
           labels = str_c(str_pad(study_long, 
                                  max_study_string+1, "right"),
                          str_pad(Chemical, 
                                  max_chem_string+1, "left"),
                          str_pad(Species, 
                                  max_species_string+3, "left")))
  
  
parameters_unicode = case_when(parameter == "mupq" ~ "expression(mu[P])",
                                          parameter == "munq" ~ "expression(mu[N])",
                                          parameter == "Knq" ~ "expression(K[N])",
                                          parameter == "psiq" ~ "expression(psi[P])",
                                          parameter == "picq" ~ "expression(pi[C])",
                                          parameter == "pimq" ~ "expression(pi[M])",
                                          parameter == "fnq" ~ "expression(f[N])",
                                          parameter == "thetaq" ~ "expression(theta)",
                                          parameter == "vq" ~ "expression(v)",
                                          parameter == "fpq" ~ "expression(f[P])")   
parameters_unicode = case_when(parameter == "mupq" ~ "\U003BC\U208P",
                                          parameter == "munq" ~ "\U003BC\U208N",
                                          parameter == "Knq" ~ "K\U208N",
                                          parameter == "psiq" ~ "\u03A8\U208P",
                                          parameter == "picq" ~ "\u03C0\U208C",
                                          parameter == "pimq" ~ "\u03C0\U208M",
                                          parameter == "fnq" ~ "f\U208N",
                                          parameter == "thetaq" ~ "\u03B8",
                                          parameter == "vq" ~ "v",
                                          parameter == "fpq" ~ "f\U208P")  
#Load all parameter sets into one dataframe
all_eec_pars <- 
  do.call(rbind,
          lapply(paste0("Agrochemical_Review/Sims/EEC/Parameter_Sets/",
                       list.files(path = "Agrochemical_Review/Sims/EEC/Parameter_Sets/")), readRDS))

#Apply r0 function to all parameter sets
  r0s <- pmap_df(all_eec_pars %>% dplyr::select(-Study_ID_Unique), r0.Ag.pars)
  
#bind parameters and results together
  all_eec_sims <- cbind(all_eec_pars, r0s) %>% 
    left_join(RFxSum)

#Separate baseline sims from all other sims
  baseline_eec_sims <- all_eec_sims %>% filter(Study_ID_Unique == "baseline_eec_parameters")
  
#Baseline summary statistics  
  base_r0_median <- median(baseline_eec_sims$R0)
  base_r0_q25 <- quantile(baseline_eec_sims$R0, 0.25)
  base_r0_q75 <- quantile(baseline_eec_sims$R0, 0.75)

  base_r0_q025 <- quantile(baseline_eec_sims$R0, 0.025)
  base_r0_q975 <- quantile(baseline_eec_sims$R0, 0.975)

#All other sims  
  eec_sims <- all_eec_sims %>% filter(Study_ID_Unique != "baseline_eec_parameters")
  
#Ggplot classic forestplot figure  #######
gg_forest_df <-  eec_sims %>% 
    group_by(Study_ID_Unique) %>% 
    summarise(med_R0 = median(R0),
              R0_25 = quantile(R0, 0.25),
              R0_75 = quantile(R0, 0.75),
              Pathway = first(Pathway),
              Class = first(Class),
              parameter = first(parameter),
              parameters_unicode = first(parameters_unicode),
              label = first(labels)) %>% 
    filter(med_R0 >= base_r0_q75 | med_R0 <= base_r0_q25) %>% 
    arrange(Class, med_R0) %>% 
    mutate(order = row_number())
  
gg_forest <- gg_forest_df %>%
    ggplot(aes(y = as.factor(order), 
               x = med_R0, xmin = R0_25, xmax = R0_75, 
               col = Pathway)) + 
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(hjust = 0, family = "mono",
                                       size = 8),
            axis.ticks.y = element_blank(),
            legend.position = "right") +
      geom_rect(xmin = base_r0_q25, xmax = base_r0_q75,
                ymin = 0, ymax = nrow(eec_sims), 
                alpha = 0.2, fill = "grey80", col = "grey80") +
      geom_vline(xintercept = base_r0_median, size = 1.2) +
      geom_point() +
      geom_errorbarh(height = 0.2) +
      scale_color_manual(values = c("bottom-up" = "#42B54099",
                                    "direct larvae" = "#925E9F99",
                                    "direct snail" = "#0099B499",
                                    "top-down" = "#ED000099")) +
      scale_y_discrete(breaks = as.factor(gg_forest_df$order), 
                       labels = gg_forest_df$label) +
      facet_grid(Class~., scales = "free", space = "free") +
      labs(x = expression(paste(R[0], " estimates", sep = " ")))
  
gg_forest

ggsave("Agrochemical_Review/Plots/Forestplot_sig_effects.pdf",
       width = 11.5, height = 8, units = "in")

#Ggplot violin forestplot figure  #######
  
  eec_sims %>% 
    ggplot(aes(y = R0, x = Study_ID_Unique, fill = Pathway)) + 
      geom_rect(ymin = base_r0_q025, ymax = base_r0_q975,
                xmin = 0, xmax = nrow(eec_sims), 
                alpha = 0.2, fill = "grey80") +
      geom_hline(yintercept = base_r0_median, size = 1.2) +
      geom_violin(alpha = 0.8) +
      scale_fill_manual(values = c("bottom-up" = "#42B54099",
                                   "direct larvae" = "#925E9F99",
                                   "direct snail" = "#0099B499",
                                   "top-down" = "#ED000099")) +
      coord_flip() +
      facet_grid(Class~., scales = "free", space = "free") +
      theme_classic()
 