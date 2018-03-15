#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(deSolve)
require(graphics)
require(tidyverse)
require(Rmisc)

#Define some plot functions ################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plot_sim_worms <- function(output.df){
  output.df %>%
    select(time = V1,  Wt =V5 , Wu = V6) %>%
    mutate(W = cov*Wt + (1-cov)*Wu) %>%
    gather(key = "Treatment", value = "Worm_Burden", -time) %>%
    ggplot(aes(x = time, y = Worm_Burden, lty = Treatment)) + theme_bw() + geom_line(col = "purple")
  
}

plot_sim_snails <- function(output.df){
  output.df %>%
    select(time = V1, S = V2, E = V3, I = V4) %>%
    gather(key = "Infection_Class", value = "Density", -time) %>%
    ggplot(aes(x = time, y = Density, col = Infection_Class)) + theme_bw() + geom_line()

}

#Load simulation arrawy containing 100 transmission parameter sets, k=0.08, single round of MDA annually #####
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/k008_mda_fullarray.Rdata")

#In full simulation array:
  #1st ([here, , , ]) dimension is time where each row represents a day
  #2nd ([ ,here, , ]) dimension is the state variables of the model including snails and treated and untreated worm burden
  #3rd ([ , ,here, ]) dimension is the simulations run over 100 candidate transmission parameter sets
  #4th ([ , , ,here]) dimension is the 4 different models considered:
    #[ , , , 1] is the simplest model with no PDD and no added NDDs (only snail carrying capacity)
    #[ , , , 2] is the model with no PDD but WITH added NDDs on worm fecundity and acquired human immunity
    #[ , , , 3] is the model WITH PDD but with no added NDDs (only snail carrying capacity)
    #[ , , , 4] is the model WITH PDD and WITH added NDDs on worm fecundity and acquired human immunity

#Some plots
plot_lowest_noPDD_noNDD <- plot_sim_worms(as.data.frame(det.runs.k008.mda1[ , , 1, 1]  ))
plot_lowest_PDD_noNDD <- plot_sim_worms(as.data.frame(det.runs.k008.mda1[ , , 1, 3]  ))
plot_lowest_noPDD_NDD <- plot_sim_worms(as.data.frame(det.runs.k008.mda1[ , , 1, 2]  ))
plot_lowest_PDD_NDD <- plot_sim_worms(as.data.frame(det.runs.k008.mda1[ , , 1, 4]  ))


#Find BBR for each model ###########
get_bbr <- function(pre_mda, post_mda){
  bbr = (1/post_mda[,c(1:(ncol(post_mda)-1))])*(pre_mda[,c(2:ncol(post_mda))] - post_mda[,c(1:(ncol(post_mda)-1))])
  mean_bbr <- colMeans(bbr)
  sd_bbr <- apply(bbr, 2, sd)
  time <- c(1:(ncol(post_mda)-1))
  
  return(t(rbind(bbr, mean_bbr, sd_bbr, time)))
}

bbr_noPDD_noNDD <- get_bbr(w.pre.k008.mda1, w.pos.k008.mda1)
bbr_PDD_noNDD <- get_bbr(w.pre.k008.mda1.pdd, w.pos.k008.mda1.pdd)
bbr_noPDD_NDD <- get_bbr(w.pre.k008.mda1.addNDD, w.pos.k008.mda1.addNDD)
bbr_PDD_NDD <- get_bbr(w.pre.k008.mda1.pdd.addNDD, w.pos.k008.mda1.pdd.addNDD)

plot_bbr <- function(bbr_object){
  as.data.frame(bbr_object) %>%
    select(mean_bbr, sd_bbr, time) %>%
    mutate(bbr_up = mean_bbr + sd_bbr,
           bbr_lo = mean_bbr - sd_bbr) %>%
    ggplot(aes(x = time, y = mean_bbr)) + geom_point() + geom_errorbar(aes(x = time, ymin = bbr_lo, ymax = bbr_up), width = 0.1) + theme_bw()
}

plot_bbr(bbr_noPDD_noNDD)
plot_bbr(bbr_PDD_noNDD)
plot_bbr(bbr_noPDD_NDD)
plot_bbr(bbr_PDD_NDD)

#Estimate slope of BBR line (aka elimination feasibility estimator) for models #######
get_bbr_slope<-function(bbr_time_series, end_time){
  if(end_time < 3) return(NA)
  else{
    coef(lm(bbr_time_series[c(1:end_time)] ~ c(1:length(bbr_time_series))[c(1:end_time)]))[2]
  }
  
}    

fill_bbr_slope <- function(bbr_object, fill_object){
  for(i in 1:nrow(bbr_object)){
    fill_object[i,] <- apply(bbr_object, 2, get_bbr_slope, end_time = i)
  }
  return(fill_object)
}

eps_noPDD_noNDD <- bbr_noPDD_noNDD
  eps_noPDD_noNDD <- fill_bbr_slope(bbr_noPDD_noNDD, eps_noPDD_noNDD)[,-c(101:103)]
  
eps_PDD_noNDD <- bbr_PDD_noNDD
  eps_PDD_noNDD <- fill_bbr_slope(bbr_PDD_noNDD, eps_PDD_noNDD)[,-c(101:103)]

eps_noPDD_NDD <- bbr_noPDD_NDD
  eps_noPDD_NDD <- fill_bbr_slope(bbr_noPDD_NDD, eps_noPDD_NDD)[,-c(101:103)]

eps_PDD_NDD <- bbr_PDD_NDD
  eps_PDD_NDD <- fill_bbr_slope(bbr_PDD_NDD, eps_PDD_NDD)[,-c(101:103)]


#Test for differences between PDD and PDD-free ###########
#include presence/absence of pdd as a regression term in the model; store p-value associated with it  

  for(m in 3:length(bbr.pdd.mean)){
    test.seq = c(1:m)
    test.df = subset(bbr.diff, time %in% test.seq)
    mod.diff = lm(bbr ~ time + pdd, data = test.df, weights = bbr.sd^-1)
    p.diff[m-2] = summary(mod.diff)$coefficients[3,4]
  }  
    
    lines(c(3:nrow(bbr.df)), p.diff, lty=3, col=3, lwd=2)

#Add final legend  
  legend('bottomright', bty='n', cex = 0.8, lty = c(1,1,2,2,3), lwd = c(2,2,1,1,3), col = c(1,2,1,2,3),
          legend = c('PDD-free model', 'PDD model', '95% CI PDD-free', '95% CI PDD', 'p-value')) 
  
#get epsilon for each sim & test if bbr slope (epsilon) significantly different from 0 over time ###########
  eps.fill = matrix(ncol = length(3:ncol(bbr)),
                    nrow = nrow(bbr))
  eps.pdd.fill = matrix(ncol = length(3:ncol(bbr)),
                        nrow = nrow(bbr))

#get epsilon (slope of BBR profile) from year 3-20 for all 100 sims    
  for(m in 1:ncol(eps.fill)){
    for(i in 1:nrow(bbr)){
      eps.fill[i,m] = lm(bbr[i,c(1:(m+2))] ~ c(1:(m+2)))$coefficients[2]
      eps.pdd.fill[i,m] = lm(bbr.pdd[i,c(1:(m+2))] ~ c(1:(m+2)))$coefficients[2]
    }
    
  } #end function
  
  #descriptive stats for epsilon estimates
    eps.mean = colMeans(eps.fill)
    eps.sd = apply(eps.fill, 2, sd)
    
    eps.pdd.mean = colMeans(eps.pdd.fill)
    eps.pdd.sd = apply(eps.pdd.fill, 2, sd)
    
#vectors to fill with p-values from t tests    
  t.fill.g0 = as.numeric()        #is pdd-free model slope greater than 0?
  t.fill.l0 = as.numeric()        #is pdd-free model slope less than 0?
  t.fill.pdd.l0 = as.numeric()    #is pdd model slope less than 0/
  t.fill.pdd.diff = as.numeric()  #are pdd and pdd-free slopes different?

  for(t in 1:ncol(eps.fill)){
    t.fill.g0[t] = t.test(x = eps.fill[,t], mu = 0, alternative = 'greater')$p.value
    t.fill.l0[t] = t.test(x = eps.fill[,t], mu = 0, alternative = 'less')$p.value
    t.fill.pdd.l0[t] = t.test(x = eps.pdd.fill[,t], mu = 0, alternative = 'less')$p.value
    t.fill.pdd.diff[t] = t.test(x = eps.fill[,t], y = eps.pdd.fill[,t], alternative = 'two.sided')$p.value
  }

#plot epsilon over time 
  
  plot(c(3:ncol(bbr)), eps.mean, pch = 16, cex = 0.8, xlim = c(0,20), ylim = c(-0.03, 0.01), xlab = 'time (yrs)',
        ylab = expression(paste('Elimination Feasibility Estimator ( ', epsilon, ')',sep = '')))
    points(c(3:ncol(bbr)), eps.pdd.mean, pch = 16, cex = 0.8, col=2)
    
    #Add error bars as st dev
    for(s in 1:length(eps.sd)){
      segments(x0 = s+2, x1 = s+2, 
               y0 = eps.mean[s] + eps.sd[s], y1 = eps.mean[s] - eps.sd[s])
      segments(x0 = s+2, x1 = s+2, 
               y0 = eps.pdd.mean[s] + eps.pdd.sd[s], y1 = eps.pdd.mean[s] - eps.pdd.sd[s], col = 2)
    }


    
#manuscript plot version ######
#worm burden trajectories
nopdd_addndd <- as.data.frame(det.runs.k008.mda1[ , , 1, 2])
nopdd_addndd <- nopdd_addndd %>%
  mutate(Model = "PDD-Free") %>%
  select(time = V1,  Wt =V5 , Wu = V6, Model = Model) %>%
  mutate(W = cov*Wt + (1-cov)*Wu,
         t_yrs = time / 365 + 60/365) %>%
  gather(key = "Treatment", value = "Worm_Burden", W)

pdd_addndd <- as.data.frame(det.runs.k008.mda1[ , , 1, 4])
  pdd_addndd <- pdd_addndd %>%
    mutate(Model = "PDD") %>%
    select(time = V1,  Wt =V5 , Wu = V6, Model = Model) %>%
    mutate(W = cov*Wt + (1-cov)*Wu,
           t_yrs = time / 365) %>%
    gather(key = "Treatment", value = "Worm_Burden", W)
            
add_ndd_comps<-rbind(pdd_addndd, nopdd_addndd)

w.ggp<- add_ndd_comps %>%
    ggplot(aes(x = t_yrs, y = Worm_Burden, color = Model)) + 
      theme_bw(base_size = 15) +
      theme(axis.title.y = element_text(size = 12)) +
      xlim(0,61) +
      ylim(0,35) +
      labs(x = 'time (yrs)',
           y = expression(paste('Mean Worm Burden (', italic('W'), ')', sep = ''))) +
      scale_color_manual(values = c('black', 'red')) +
      geom_line() +
      annotate('text', x = 0.25, y = 35, label = 'a', size = 5)

#Bounce back rates  
bbr_PDD_NDD_gg <- as.data.frame(bbr_PDD_NDD) %>%
  select(101:103) %>%
  mutate(Model = "PDD")

bbr_noPDD_NDD_gg <- as.data.frame(bbr_noPDD_NDD) %>%
  select(101:103) %>%
  mutate(Model = "PDD-Free")

bbr.gg = rbind(bbr_PDD_NDD_gg, bbr_noPDD_NDD_gg)

bbr.ggp <- bbr.gg %>%
  ggplot(aes(x = time, y = mean_bbr, color = Model)) +
    theme_bw(base_size = 15) +
    theme(legend.position = 'none', axis.title.y = element_text(size = 12)) +
    xlim(0,20) +
    labs(x = 'time (yrs)',
         y = expression(paste('Bounce Back Rate (', italic('BBR'), ')', sep = ''))) +
    scale_color_manual(values = alpha(c( 'black', 'red'), .9)) +
    geom_point() +
    geom_errorbar(aes(x = time, ymin = (mean_bbr - sd_bbr), ymax = (mean_bbr + sd_bbr)), 
                  width = 0.1) +
    annotate('text', x = 0, y = 0.34, label = 'b', size = 5)

#epsilon estimates
eps_noPDD_NDD <- as.data.frame(eps_noPDD_NDD) %>%
  mutate(mean_eps = rowMeans(.),
         sd_eps = apply(., 1, sd),
         time = c(1:19),
         Model = "PDD-Free") %>%
  select(mean_eps, sd_eps, time, Model)

eps_PDD_NDD <- as.data.frame(eps_PDD_NDD) %>%
  mutate(mean_eps = rowMeans(.),
         sd_eps = apply(., 1, sd),
         time = c(1:19),
         Model = "PDD") %>%
  select(mean_eps, sd_eps, time, Model)

eps.gg <- rbind(eps_noPDD_NDD, eps_PDD_NDD)

eps.ggp <- eps.gg %>%
  ggplot(aes(x = time, y = mean_eps, fill = Model)) +
    theme_bw(base_size = 15) +
    theme(legend.position = 'none', axis.title.y = element_text(size = 10)) +
    xlim(0,20) +
    ylim(-0.035,0.01) +
    labs(x = 'time (yrs)',
         y = expression(paste('Elimination Feasibility Coefficient (',epsilon, ')', sep = ''))) +
    geom_bar(position = 'dodge', stat = 'identity', width = 0.8) +
    scale_fill_manual(values = alpha(c('black', 'red'), .6)) +
    geom_errorbar(aes(x = time, ymin = (mean_eps - sd_eps), ymax = (mean_eps + sd_eps)), 
                  width = 0.1, position = position_dodge(width = 0.8)) +
    geom_vline(xintercept = 4.5, lty=2) +
    annotate('text', x = 0, y = 0.009, label = 'c', size = 5)


#combine plots with multiplot
windows(width = 14, height = 9)
  multiplot(w.ggp, bbr.ggp, eps.ggp, layout = matrix(c(1,1,2,3), nrow = 2, byrow = T))
graphics.off()  
  
#Save as tiff
tiff("Elimination_Feasibility/plots/PLoS_Figs_PostReview/Fig3.tiff", height = 14, width = 19.05, units = 'cm', 
     compression = "lzw", res = 300) 
  multiplot(w.ggp, bbr.ggp, eps.ggp, layout = matrix(c(1,1,2,3), nrow = 2, byrow = T))
  
dev.off()
  
