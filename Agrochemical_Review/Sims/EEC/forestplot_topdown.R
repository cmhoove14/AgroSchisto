#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

library(tidyverse)
library(pals)
library(forestplot)

#GGplot theme for manuscripts
source("Agrochemical_Review/Sims/ggplot_theme.R")

#load R0 function
source("Agrochemical_Review/Models/r0_of_q.R")

#Response functions summary
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

rfx_dir <- "Agrochemical_Review/Response_Fxs/"

rfx_files <- list.files(path = "Agrochemical_Review/Response_Fxs",
                        pattern = "_fit.R", recursive = TRUE)

#Initially was trying to work on full data frame, but introduces issues since there are so many objects required, therefore separating jobs by pathway. Also makse sense since we're constructing separate forestplots for each pathway
rfx_topdown <- rfx_sum %>% filter(Pathway == "top-down") %>% rename(Parameter = parameter,
                                                                    Study = study_long,
                                                                    study_abrev = Study)

#Load all top-down response functions ######
topdown_studies <- unique(rfx_topdown$study_abrev)

sapply(paste0(rfx_dir, rfx_files[unlist(sapply(topdown_studies, grep, rfx_files, ignore.case = TRUE))]), source)


#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for muPq response functions
muPq_r0 <- function(rfx, eec){
  print(rfx)
  muPq_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     #muPqx = 0,
                     fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(muPq_r0s),
                    r0_025 = quantile(muPq_r0s, 0.25),
                    r0_975 = quantile(muPq_r0s, 0.75)))
}

rfx_topdown_muPq <- rfx_topdown %>% filter(Parameter == "muPq" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, muPq_r0))

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for psiq response functions
psiq_r0 <- function(rfx, eec){
  print(rfx)
  psiq_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     muPqx = 0,
                     fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     #psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(psiq_r0s),
                    r0_025 = quantile(psiq_r0s, 0.25),
                    r0_975 = quantile(psiq_r0s, 0.75)))
}

rfx_topdown_psiq <- rfx_topdown %>% filter(Parameter == "psiq" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, psiq_r0))

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for fPq response functions
fPq_r0 <- function(rfx, eec){
  print(rfx)
  fPq_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     muPqx = 0,
                     #fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(fPq_r0s),
                    r0_025 = quantile(fPq_r0s, 0.25),
                    r0_975 = quantile(fPq_r0s, 0.75)))
}

rfx_topdown_fPq <- rfx_topdown %>% filter(Parameter == "fPq" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, psiq_r0))

#Put dfs back together
rfx_topdown_all <- rbind(rfx_topdown_muPq, rfx_topdown_psiq, rfx_topdown_fPq) %>% 
  mutate(r0_med_rel = (r0_med / r0.fix()[3]) * 100,
         r0_025_rel = (r0_025 / r0.fix()[3]) * 100,
         r0_975_rel = (r0_975 / r0.fix()[3]) * 100)

save(rfx_topdown_all, file = "Agrochemical_Review/Sims/EEC/all_topdown_summary.RData")

#make forestplot #######
fp_studies <- c("Study",rfx_topdown_all %>% arrange(Parameter, Chemical, Study) %>% pull(Study))
fp_chems <- c("Agrochemical", rfx_topdown_all %>% arrange(Parameter, Chemical, Study) %>% pull(Chemical))
fp_species <- c("Species",rfx_topdown_all %>% arrange(Parameter, Chemical, Study) %>% pull(Species))
fp_pars <- c("Parameter",rfx_topdown_all %>% arrange(Parameter, Chemical, Study) %>% pull(Parameter))

fp_text <- list(as.list(fp_studies),
                as.list(fp_chems),
                as.list(fp_species),
                as.list(fp_pars))

  forestplot(labeltext = fp_text, 
             #fn.ci_norm = c(fpDrawNormalCI),
             line.margin = 0.1,
             mean = matrix(c(NA, rfx_topdown_all %>% arrange(Parameter, Chemical, Study) %>% pull(r0_med))),
             lower =matrix(c(NA, rfx_topdown_all %>% arrange(Parameter, Chemical, Study) %>% pull(r0_025))),
             upper =matrix(c(NA, rfx_topdown_all %>% arrange(Parameter, Chemical, Study) %>% pull(r0_975))),
             new_page = TRUE,
             is.summary = c(TRUE, rep(FALSE, nrow(rfx_topdown_all))),
             hrzl_lines = list('2' = gpar(col = 'grey')),
             txt_gp = fpTxtGp(xlab = gpar(cex = 1.2),
                              ticks = gpar(cex = 1.1)),
             vertices = TRUE,
             zero = r0.fix()[3],
             boxsize = 0.4,
             col = fpColors(box = 'black', lines = 'black'),
             clip=c(-Inf,Inf),
             xlab = expression(paste(R['0'])),
             xticks = c(0, 1,2,3,4))


#Forestplot-ish thing with ggplot ###########
my_labs <- list(bquote(mu[P]),bquote(Psi),bquote(f[P]))

rfx_topdown_all %>% 
  ggplot(aes(x = reorder(Chemical, desc(Chemical)), y = r0_med, shape = Parameter, col = Study)) + 
    geom_hline(yintercept = r0.fix()[3], lty = 2) +
    geom_point(size = 2.5, position = position_dodge(0.9)) + 
    geom_errorbar(aes(ymin = r0_025, ymax = r0_975, x = Chemical, width = 0.01), 
                  position = position_dodge(0.9)) +
    theme_ms() + ylim(0,4) + coord_flip() + ylab(expression(paste(R['0']))) +
    ggtitle("Top-down effects") + xlab("Chemical") +
    scale_color_manual(values = glasbey()) + 
    scale_shape_manual(values = c(17,15,16), breaks = c("muPq", "psiq", "fPq"),
                       labels = my_labs) 

ggsave(paste('~/RemaisWork/Schisto/Agro_Review/Figures/EEC_forest/ggplot_forest_topdown', Sys.Date(), '.tiff', sep = ''),
        width = 7.3, height = 7.3, dpi = 600)
