require(tidyverse)
require(pals)

#GGplot theme for manuscripts
source("Agrochemical_Review/Sims/ggplot_theme.R")

#r0 function
source("Agrochemical_Review/Models/r0_of_q.R")

#Load all response function summaries
load("Agrochemical_Review/Sims/EEC/all_topdown_summary.RData")
load("Agrochemical_Review/Sims/EEC/all_dirsnail_summary.RData")
load("Agrochemical_Review/Sims/EEC/all_dirlarv_summary.RData")
load("Agrochemical_Review/Sims/EEC/all_bottomup_summary.RData")

rfx_all_sum <- rbind(rfx_bottomup_all, rfx_dirlarv_all, rfx_dirsnail_all, rfx_topdown_all)

#Forestplot with all insecticide-driven effects
#parameter labels
ins_labs <- list(bquote(paste("Snail reproduction (", f[N], ")")),
                 bquote(paste("Predator reproduction (", f[P], ")")),
                  bquote(paste("Snail mortality (", mu[N], ")")), 
                  bquote(paste("Predator mortality (", mu[P], ")")),
                  bquote(paste("Cercarial mortality (", pi[C], ")")), 
                  bquote(paste("Miracidial mortality (", pi[M], ")")), 
                  bquote(paste("Predator attack rate (", psi, ")")), 
                  bquote(paste("Egg viability (", v, ")")))

ins_forest <- rfx_all_sum %>% 
  filter(System != "Other" & Class == "Insecticide") %>% 
    ggplot(aes(x = reorder(Study, desc(Study)), y = r0_med, col = Chemical, shape = Parameter)) + 
      geom_hline(yintercept = r0.fix()[3], lty = 2) +
      geom_point(size = 3, position = position_dodge(0.9)) + 
      geom_errorbar(aes(ymin = r0_025, ymax = r0_975, x = Study, col = Chemical, width = 0.01), 
                    position = position_dodge(0.9)) +
      theme_ms() + ylim(0,4) + coord_flip() + ylab(expression(paste(R['0']))) + xlab("") +
      ggtitle("Summary of insecticide effects on components of schistosomiasis transmission") + 
      scale_color_manual(values = glasbey()) + 
      scale_shape_manual(values = c(0, 1, 15, 16, 
                                    2, 17, 
                                    12, 4), 
                         breaks = c("fNq", "fPq","muNq", "muPq",
                                    "piC", "piM", 
                                    "psiq", "vq"),
                         labels = ins_labs)

ggsave(paste('~/RemaisWork/Schisto/Agro_Review/Figures/EEC_forest/ggplot_forest_all_ins_', Sys.Date(), '.png', sep = ''),
       width = 11, height = 8.5, dpi = 600)

#Forestplot with all herbicide and fertilizer-driven effects
#parameter labels
herb_labs <- list(bquote(paste("Snail reproduction (", f[N], ")")),
                  bquote(paste("Snail mortality (", mu[N], ")")), 
                  bquote(paste("Predator mortality (", mu[P], ")")),
                  bquote(paste("Snail carrying capacity (", phi[N], ")")),
                  bquote(paste("Cercarial mortality (", pi[C], ")")), 
                  bquote(paste("Miracidial mortality (", pi[M], ")")), 
                  bquote(paste("Predator attack rate (", psi, ")")), 
                  bquote(paste("Egg viability (", v, ")")))

rfx_all_sum %>% 
  filter(System != "Other" & Class %in% c("Herbicide", "Fertilizer")) %>% 
    ggplot(aes(x = reorder(Study, desc(Study)), y = r0_med, col = Chemical, shape = Parameter)) + 
      geom_hline(yintercept = r0.fix()[3], lty = 2) +
      geom_point(size = 3, position = position_dodge(0.9)) + 
      geom_errorbar(aes(ymin = r0_025, ymax = r0_975, x = Study, col = Chemical, width = 0.01), 
                    position = position_dodge(0.9)) +
      theme_ms() + ylim(0,4) + coord_flip() + ylab(expression(paste(R['0']))) + xlab("") +
      ggtitle("Summary of herbicide and fertilizer effects on components of schistosomiasis transmission") + 
      scale_color_manual(values = glasbey()) + 
      scale_shape_manual(values = c(0, 15, 16, 
                                    13, 2, 17, 
                                    12, 4), 
                         breaks = c("fNq", "muNq", "muPq",
                                    "phiNq", "piC", "piM", 
                                    "psiq", "vq"),
                         labels = herb_labs)

ggsave(paste('~/RemaisWork/Schisto/Agro_Review/Figures/EEC_forest/ggplot_forest_all_herbs&ferts_', Sys.Date(), '.png', sep = ''),
       width = 11, height = 8.5, dpi = 600)
