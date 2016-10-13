require(ggplot2)
require(drc)

#Atrazine effect on phi_Nq (snail carrying capacity) ################
atra.df<-data.frame('atra' = c(0,1,10,30,100),              #Raw atrazine concentration (ppb)
                    'logatra' = log(c(0,1,10,30,100)+1),    #Log atrazine concentration (ppb)
                    'phiNq' = c(0,0.2888,0.6535,0,1.3215))  #Snail population response measured as peak snail growth rate and 
                                                            #interpreted as changes in snail carrying capacity

plot(atra.df$logatra, atra.df$phiNq, pch = 16)

atra_mod<-glm(phiNq ~ logatra+0, data=atra.df)

atra.slope = atra_mod$coefficients[1]

#24-hr mortality functions for predator (macrobrachium) populations from Satapornvanit et al chemosphere paper
  sap.mort<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_mort.csv")
    ggplot(data = sap.mort, aes(x = log(conc+1), y = mort, color = chem)) +
        theme_bw() +
        labs(y = 'Mortality (%)', x = 'log+1 Chem Conc (ppb)', title = 'M. rosenbergii 24 hr mortality') +
        xlim(0,12.5) +
        ylim(0,100) +
        geom_point(size = 4) 
    
    z.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'zinc']/100 ~ sap.mort$conc[sap.mort$chem == 'zinc'],
                    data = sap.mort, type = 'binomial', fct = LL.2())
    
    ch.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'chlorpyrifos']/100 ~ sap.mort$conc[sap.mort$chem == 'chlorpyrifos'],
                     data = sap.mort, type = 'binomial', fct = LL.2())
    
    dim.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'dimethoate']/100 ~ sap.mort$conc[sap.mort$chem == 'dimethoate'],
                      data = sap.mort, type = 'binomial', fct = LL.2())
    
    pr.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'profenofos']/100 ~ sap.mort$conc[sap.mort$chem == 'profenofos'],
                     data = sap.mort, type = 'binomial', fct = LL.2())
    
    #NOTE: carbendazim not modeled because it has no effect on toxicity, but will be included as such in simulations
    
  sap.fr<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_feed_rate.csv")
    fr.z = subset(sap.fr, chem == 'zinc' & feed_rate > 0 )
    fr.ch = subset(sap.fr, chem == 'chlorpyrifos' & feed_rate > 0)

    #Get slope only and fix y intercept at 1 (i.e. no reduction in feeding rate)
      z.mod.fr<-lm(log(per_control) ~ conc + 0, data = fr.z) 
      ch.mod.fr<-lm(log(per_control) ~ conc + 0, data = fr.ch) 

  ggplot(data = fr.z, aes(x = conc, y = log(per_control))) +
      theme_bw() +
      labs(y = 'Log Feed rate (%control)', x = 'Zinc Conc (ppb)', title = 'M. rosenbergii feed reduction') +
      xlim(0,max(fr.z$conc)) +
      ylim(-3,0) +
      geom_abline(intercept = 0, slope = z.mod.fr$coefficients, lty = 2, color = 'red') +
      geom_point(size = 4) +
      geom_errorbar(aes(x = conc, ymax = log(per_max), ymin = log(per_min)), width = 0.1)
    
  ggplot(data = fr.ch, aes(x = conc, y = log(per_control))) +
    theme_bw() +
    labs(y = 'Log Feed rate (%control)', x = 'Chlorpyrifos Conc (ppb)', title = 'M. rosenbergii feed reduction') +
    xlim(0,max(fr.ch$conc)) +
    ylim(-3,0) +
    geom_abline(intercept = 0, slope = ch.mod.fr$coefficients, lty = 2, color = 'red') +
    geom_point(size = 4) +
    geom_errorbar(aes(x = conc, ymax = log(per_max), ymin = log(per_min)), width = 0.01)