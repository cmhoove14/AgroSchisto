#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


require(ggplot2)
require(drc)
require(LW1949)

#Insecticide toxicity to crustaceans from Halstead et al; 4-day mortality endpoints ################
  data<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Cray.LC50.2012.csv')
  
#Malathion ******************************************************************************########
  mal<-subset(data, Chem == 'Mal')
    mal.sum = data.frame(mal.c = unique(mal$Conc),
                         mal.total = 5,
                         mal.d = 0,
                         mort = 0)
  #DRC package analysis  
    for(i in 1:length(mal.sum$mal.d)){
      mal.sum$mal.d[i] = sum(mal$Dead[mal$Conc == mal.sum$mal.c[i]])
    }
    
    mal.sum$mort = mal.sum$mal.d / mal.sum$mal.total
    
      mal.mod<-drm(mal.d / mal.total ~ mal.c, weights = mal.total, data = mal.sum, type = 'binomial',  
                   fct = LL.2(names = c("b", "e"),
                              fixed = c(NA, NA)))

     muPq_mal_Halstead<-function(In){
        predict(mal.mod, data.frame(mal.c = In), interval = 'confidence', level = 0.95)
      }  
     
     plot(mal.sum$mal.c, mal.sum$mort, pch = 16, ylim = c(0,1),
          xlab = 'Malathion', ylab = 'Mortality', 
          main = expression(paste('Malathion toxicity to ', italic('P. clarkii')))) 
     
      lines(seq(0, max(mal$Conc), 100), sapply(seq(0, max(mal$Conc), 100), muPq_mal_Halstead, simplify = T)[1,],
            lty = 2, col = 2)
      lines(seq(0, max(mal$Conc), 100), sapply(seq(0, max(mal$Conc), 100), muPq_mal_Halstead, simplify = T)[2,],
            lty = 3, col = 2)
      lines(seq(0, max(mal$Conc), 100), sapply(seq(0, max(mal$Conc), 100), muPq_mal_Halstead, simplify = T)[3,],
            lty = 3, col = 2)
   
    par.mal = c(coef(mal.mod), 'Lower Limit:(Intercept)' = 0, 'Upper Limit:(Intercept)' = 1)[c(1,3,4,2)]  
           
      muPq_mal_Halstead_uncertainty<-function(In){
        rdrm(1, LL.2(), coef(mal.mod), In, yerror = 'rbinom', ypar = 5)$y / 5 #estimate deaths / live 
      }
        points(seq(0, 40000, 40), sapply(seq(0, 40000, 40), muPq_mal_Halstead_uncertainty),
             pch = 5, col=4, cex = 0.5)
  #LW1949 analysis
    mal.lw1949 = dataprep(dose = mal.sum$mal.c, ntot = mal.sum$mal.total, nfx = mal.sum$mal.d)
    slp.mal = fitLWauto(mal.lw1949)
    LW.mal = LWestimate(slp.mal, mal.lw1949)
  #Doesn't work because of lack of doses that have an effect    
    
#Chlorpyrifos   *************************************************************************############
  chlor<-subset(data, Chem == 'Chlor')
    chlor.sum = data.frame(chlor.c = unique(chlor$Conc),
                         chlor.total = 5,
                         chlor.d = 0,
                         mort = 0)
      
    for(i in 1:length(chlor.sum$chlor.d)){
        chlor.sum$chlor.d[i] = sum(chlor$Dead[chlor$Conc == chlor.sum$chlor.c[i]])
    }
      
    chlor.sum$mort = chlor.sum$chlor.d / chlor.sum$chlor.total
      
    chlor.mod<-drm(chlor.d / chlor.total ~ chlor.c, weights = chlor.total, data = chlor.sum, type = 'binomial',  
                   fct = LL.2(names = c('b', 'e'),
                              fixed = c(NA, NA)))

    muPq_chlor_Halstead<-function(In){
      predict(chlor.mod, data.frame(chlor.c = In), interval = 'confidence', level = 0.95)
    }  
      
    plot(chlor.sum$chlor.c, chlor.sum$mort, pch = 16, ylim = c(0,1),
         xlab = 'Chlorpyrifos', ylab = 'Mortality', 
         main = expression(paste('Chlorpyrifos toxicity to ', italic('P. clarkii')))) 
    
    lines(seq(0, max(chlor$Conc), 1), sapply(seq(0, max(chlor$Conc), 1), muPq_chlor_Halstead, simplify = T)[1,],
          lty = 2, col = 2)
    lines(seq(0, max(chlor$Conc), 1), sapply(seq(0, max(chlor$Conc), 1), muPq_chlor_Halstead, simplify = T)[2,],
          lty = 3, col = 2)
    lines(seq(0, max(chlor$Conc), 1), sapply(seq(0, max(chlor$Conc), 1), muPq_chlor_Halstead, simplify = T)[3,],
          lty = 3, col = 2)
     
      muPq_chlor_Halstead_uncertainty_DRC<-function(In){
        rdrm(1, LL.2(), coef(chlor.mod), In, yerror = 'rbinom', ypar = 5)$y / 5
      } 
        
      points(c(0:100), sapply(c(0:100), muPq_chlor_Halstead_uncertainty_DRC),
           pch = 5, col=4, cex = 0.5)
      
    #LW1949 analysis
      chlor.lw1949 = dataprep(dose = chlor.sum$chlor.c, ntot = chlor.sum$chlor.total, nfx = chlor.sum$chlor.d)
        chlor.lw1949$cbitpfx <- constrain(chlor.lw1949$bitpfx, probit(c(5e-04, 0.9995)))
        chlor.lwmod = lm(cbitpfx ~ log10dose, data = chlor.lw1949[chlor.lw1949$LWkeep, ])
        
      lines(seq(0, max(chlor$Conc), 1), 
            pnorm(predict(chlor.lwmod, newdata = data.frame(log10dose = log10(seq(0, max(chlor$Conc), 1))), 
                          interval = 'confidence', level = 0.95)[,1]), lty = 2, col = 4)  

      lines(seq(0, max(chlor$Conc), 1), 
            pnorm(predict(chlor.lwmod, newdata = data.frame(log10dose = log10(seq(0, max(chlor$Conc), 1))), 
                          interval = 'confidence', level = 0.95)[,2]), lty = 3, col = 4) 
      lines(seq(0, max(chlor$Conc), 1), 
            pnorm(predict(chlor.lwmod, newdata = data.frame(log10dose = log10(seq(0, max(chlor$Conc), 1))), 
                          interval = 'confidence', level = 0.95)[,3]), lty = 3, col = 4) 
      
      muPq_chlor_Halstead_uncertainty = function(In){
        if(In==0) mup = 0 else{
          init = predict(chlor.lwmod, newdata = data.frame(log10dose = log10(In)), se.fit = T)
          mup = pnorm(rnorm(1, init$fit, init$se.fit))#
        }
        return(mup)
      }
      
      points(seq(0,70,0.5), sapply(seq(0,70,0.5), muPq_chlor_Halstead_uncertainty),
             pch = 5, col=6, cex = 0.5)

#Terbufos   *****************************************************************************########
  terb<-subset(data, Chem == 'Terb')
    terb.sum = data.frame(terb.c = unique(terb$Conc),
                             terb.total = 5,
                             terb.d = 0,
                             mort = 0)
      
      for(i in 1:length(terb.sum$terb.d)){
        terb.sum$terb.d[i] = sum(terb$Dead[terb$Conc == terb.sum$terb.c[i]])
      }
      
      terb.sum$mort = terb.sum$terb.d / terb.sum$terb.total
      
  terb.mod<-drm(terb.d / terb.total ~ terb.c, weights = terb.total, data = terb.sum, type = 'binomial', 
                fct = LL.2(names = c('b', 'e'),
                           fixed = c(NA, NA)))
  
  muPq_terb_Halstead<-function(In){
    predict(terb.mod, newdata = data.frame(terb.c = In), interval = 'confidence', level = 0.95)
  }
        
    plot(terb.sum$terb.c, terb.sum$mort, pch = 16, xlab = 'Terbufos Concentration (ppb)', 
         ylab = 'prop dead', ylim = c(0,1),
         main = expression(paste('Terbufos toxicity to ', italic('P. clarkii'))))
      lines(c(0:170), sapply(c(0:170), muPq_terb_Halstead)[1,], col = 2, lty=2)
      lines(c(0:170), sapply(c(0:170), muPq_terb_Halstead)[2,], col = 2, lty=3)
      lines(c(0:170), sapply(c(0:170), muPq_terb_Halstead)[3,], col = 2, lty=3)
      

    muPq_terb_Halstead_uncertainty_DRC<-function(In){
      rdrm(1, LL.2(), coef(terb.mod), In, yerror = 'rbinom', ypar = 5)$y / 5
    } 
    
    points(c(0:170), sapply(c(0:170), muPq_terb_Halstead_uncertainty_DRC),
           pch = 5, col=4, cex = 0.5)

#LW1949 analysis
terb.lw1949 = dataprep(dose = terb.sum$terb.c, ntot = terb.sum$terb.total, nfx = terb.sum$terb.d)
  terb.lw1949$cbitpfx <- constrain(terb.lw1949$bitpfx, probit(c(5e-04, 0.9995)))
  terb.lwmod = lm(cbitpfx ~ log10dose, data = terb.lw1949[terb.lw1949$LWkeep, ])

  lines(seq(0, max(terb$Conc), 1), 
        pnorm(predict(terb.lwmod, newdata = data.frame(log10dose = log10(seq(0, max(terb$Conc), 1))), 
                      interval = 'confidence', level = 0.95)[,1]), lty = 2, col = 4)  
  
  lines(seq(0, max(terb$Conc), 1), 
        pnorm(predict(terb.lwmod, newdata = data.frame(log10dose = log10(seq(0, max(terb$Conc), 1))), 
                      interval = 'confidence', level = 0.95)[,2]), lty = 3, col = 4) 
  lines(seq(0, max(terb$Conc), 1), 
        pnorm(predict(terb.lwmod, newdata = data.frame(log10dose = log10(seq(0, max(terb$Conc), 1))), 
                      interval = 'confidence', level = 0.95)[,3]), lty = 3, col = 4) 

muPq_terb_Halstead_uncertainty = function(In){
  if(In==0) mup = 0 else{
    init = predict(terb.lwmod, newdata = data.frame(log10dose = log10(In)), se.fit = T)
    mup = pnorm(rnorm(1, init$fit, init$se.fit))# daily
  }
  return(mup)
}

points(seq(0,170,0.5), sapply(seq(0,170,0.5), muPq_terb_Halstead_uncertainty),
       pch = 5, col=6, cex = 0.5)

        
#Lambda-cyhalothrin   *******************************************************************##########
  lamcy<-subset(data, Chem == 'Lambda')
  lamcy.sum = data.frame(lamcy.c = unique(lamcy$Conc),
                        lamcy.total = 5,
                        lamcy.d = 0,
                        mort = 0)
  
    for(i in 1:length(lamcy.sum$lamcy.d)){
      lamcy.sum$lamcy.d[i] = sum(lamcy$Dead[lamcy$Conc == lamcy.sum$lamcy.c[i]])
    }
    
    lamcy.sum$mort = lamcy.sum$lamcy.d / lamcy.sum$lamcy.total
  
  lamcy.mod<-drm(lamcy.d / lamcy.total ~ lamcy.c, weights = lamcy.total, data = lamcy.sum, type = 'binomial', 
                fct = LL.2(names = c('b', 'e'),
                           fixed = c(NA, NA)))
 
    muPq_lamcy_Halstead<-function(In){
      predict(lamcy.mod, data.frame(lamcy.c = In), interval = 'confidence', level = 0.95)
    }  
  
  
  plot(lamcy.sum$lamcy.c, lamcy.sum$mort, pch = 16, xlab = expression(paste(lambda, '-cyhalothrin (ppb)')), 
       ylab = expression(paste(mu[P])), ylim = c(0,1),
       main = expression(paste(lambda,'-cyhalothrin toxicity to ', italic('P. clarkii'))))
    lines(seq(0,3,0.01), sapply(seq(0,3,0.01), muPq_lamcy_Halstead)[1,], col = 2, lty=2)
    lines(seq(0,3,0.01), sapply(seq(0,3,0.01), muPq_lamcy_Halstead)[2,], col = 2, lty=3)
    lines(seq(0,3,0.01), sapply(seq(0,3,0.01), muPq_lamcy_Halstead)[3,], col = 2, lty=3)
  
  muPq_lamcy_Halstead_uncertainty_DRC<-function(In){
    rdrm(1, LL.2(), coef(lamcy.mod), In, yerror = 'rbinom', ypar = 5)$y / 5
  } 
  
  points(seq(0, 3, 0.03), sapply(seq(0, 3, 0.03), muPq_lamcy_Halstead_uncertainty_DRC),
       pch = 5, col=4, cex = 0.5)

#LW1949 analysis
lamcy.lw1949 = dataprep(dose = lamcy.sum$lamcy.c, ntot = lamcy.sum$lamcy.total, nfx = lamcy.sum$lamcy.d)
lamcy.lw1949$cbitpfx <- constrain(lamcy.lw1949$bitpfx, probit(c(5e-04, 0.9995)))
lamcy.lwmod = lm(cbitpfx ~ log10dose, data = lamcy.lw1949[lamcy.lw1949$LWkeep, ])

  lines(seq(0,3,0.003), 
        pnorm(predict(lamcy.lwmod, newdata = data.frame(log10dose = log10(seq(0,3,0.003))), 
                      interval = 'confidence', level = 0.95)[,1]), lty = 2, col = 4)  
  
  lines(seq(0,3,0.003), 
        pnorm(predict(lamcy.lwmod, newdata = data.frame(log10dose = log10(seq(0,3,0.003))), 
                      interval = 'confidence', level = 0.95)[,2]), lty = 3, col = 4) 
  lines(seq(0,3,0.003), 
        pnorm(predict(lamcy.lwmod, newdata = data.frame(log10dose = log10(seq(0,3,0.003))), 
                      interval = 'confidence', level = 0.95)[,3]), lty = 3, col = 4) 

  muPq_lamcy_Halstead_uncertainty = function(In){
    if(In==0) mup = 0 else{
      init = predict(lamcy.lwmod, newdata = data.frame(log10dose = log10(In)), se.fit = T)
      mup = pnorm(rnorm(1, init$fit, init$se.fit))# daily
    }
    return(mup)
  }
  
    points(seq(0,3,0.03), sapply(seq(0,3,0.03), muPq_lamcy_Halstead_uncertainty),
           pch = 5, col=6, cex = 0.5)
#esfenvalerate  *************************************************************************#######
  esfen<-subset(data, Chem == 'Esfen')
    esfen.sum = data.frame(esfen.c = unique(esfen$Conc),
                           esfen.total = 5,
                           esfen.d = 0,
                           mort = 0)
  
      for(i in 1:length(esfen.sum$esfen.d)){
        esfen.sum$esfen.d[i] = sum(esfen$Dead[esfen$Conc == esfen.sum$esfen.c[i]])
      }
      
      esfen.sum$mort = esfen.sum$esfen.d / esfen.sum$esfen.total
  
  esfen.mod<-drm(esfen.d / esfen.total ~ esfen.c, weights = esfen.total, data = esfen.sum, type = 'binomial', 
                 fct = LL.2(names = c('b', 'e'),
                            fixed = c(NA, NA)))

    muPq_esfen_Halstead<-function(In){
      predict(esfen.mod, data.frame(esfen.c = In), interval = 'confidence', level = 0.95)
    }  
  
  plot(esfen.sum$esfen.c, esfen.sum$mort, pch = 16, xlab = 'Esfenvalerate Concentration (ppb)', 
       ylab = expression(paste(mu[P])), ylim = c(0,1),
       main = expression(paste('Esfenvalerate toxicity to ', italic('P. clarkii'))))
    lines(seq(0,25,0.1), sapply(seq(0,25,0.1), muPq_esfen_Halstead)[1,], col = 2, lty=2)
    lines(seq(0,25,0.1), sapply(seq(0,25,0.1), muPq_esfen_Halstead)[2,], col = 2, lty=3)
    lines(seq(0,25,0.1), sapply(seq(0,25,0.1), muPq_esfen_Halstead)[3,], col = 2, lty=3)
    
  muPq_esfen_Halstead_uncertainty_DRC<-function(In){
    rdrm(1, LL.2(), coef(esfen.mod), In, yerror = 'rbinom', ypar = 5)$y / 5 
  } 
  
  points(seq(0, 25, 0.1), sapply(seq(0, 25, 0.1), muPq_esfen_Halstead_uncertainty_DRC),
         pch = 5, col=4, cex = 0.5)

#LW1949 analysis
esfen.lw1949 = dataprep(dose = esfen.sum$esfen.c, ntot = esfen.sum$esfen.total, nfx = esfen.sum$esfen.d)
esfen.lw1949$cbitpfx <- constrain(esfen.lw1949$bitpfx, probit(c(5e-04, 0.9995)))
esfen.lwmod = lm(cbitpfx ~ log10dose, data = esfen.lw1949[esfen.lw1949$LWkeep, ])

  lines(seq(0, 25, 0.1), 
        pnorm(predict(esfen.lwmod, newdata = data.frame(log10dose = log10(seq(0, 25, 0.1))), 
                      interval = 'confidence', level = 0.95)[,1]), lty = 2, col = 4)  
  
  lines(seq(0, 25, 0.1), 
        pnorm(predict(esfen.lwmod, newdata = data.frame(log10dose = log10(seq(0, 25, 0.1))), 
                      interval = 'confidence', level = 0.95)[,2]), lty = 3, col = 4) 
  lines(seq(0, 25, 0.1), 
        pnorm(predict(esfen.lwmod, newdata = data.frame(log10dose = log10(seq(0, 25, 0.1))), 
                      interval = 'confidence', level = 0.95)[,3]), lty = 3, col = 4) 

  muPq_esfen_Halstead_uncertainty = function(In){
    if(In==0) mup = 0 else{
      init = predict(esfen.lwmod, newdata = data.frame(log10dose = log10(In)), se.fit = T)
      mup = pnorm(rnorm(1, init$fit, init$se.fit))# daily
    }
    return(mup)
  }
  
    points(seq(0, 25, 0.1), sapply(seq(0, 25, 0.1), muPq_esfen_Halstead_uncertainty),
           pch = 5, col=6, cex = 0.5)

#Permethrin *****************************************************************************##########
  perm<-subset(data, Chem == 'Perm')
    perm.sum = data.frame(perm.c = unique(perm$Conc),
                           perm.total = 5,
                           perm.d = 0,
                           mort = 0)
  
      for(i in 1:length(perm.sum$perm.d)){
        perm.sum$perm.d[i] = sum(perm$Dead[perm$Conc == perm.sum$perm.c[i]])
      }
      
      perm.sum$mort = perm.sum$perm.d / perm.sum$perm.total
  
  perm.mod<-drm(perm.d / perm.total ~ perm.c, weights = perm.total, data = perm.sum, type = 'binomial', 
                 fct = LL.2(names = c('b', 'e'),
                            fixed = c(NA, NA)))
  
  muPq_perm_Halstead<-function(In){
    predict(perm.mod, data.frame(perm.c = In), interval = 'confidence', level = 0.95)
  }  
  
  plot(perm.sum$perm.c, perm.sum$mort, pch = 16, xlab = 'Permethrin Concentration (ppb)', 
       ylab = expression(paste(mu[P])), ylim = c(0,1),
       main = expression(paste('Permethrin toxicity to ', italic('P. clarkii'))))
    lines(seq(0,6,0.1), sapply(seq(0,6,0.1), muPq_perm_Halstead)[1,], col = 2, lty=2)
    lines(seq(0,6,0.1), sapply(seq(0,6,0.1), muPq_perm_Halstead)[2,], col = 2, lty=3)
    lines(seq(0,6,0.1), sapply(seq(0,6,0.1), muPq_perm_Halstead)[3,], col = 2, lty=3)
    
  par.perm = c(coef(perm.mod), 'Lower Limit:(Intercept)' = 0, 'Upper Limit:(Intercept)' = 1)[c(1,3,4,2)]
  
  muPq_perm_Halstead_uncertainty_DRC<-function(In){
    rdrm(1, LL.2(), coef(perm.mod), In, yerror = 'rbinom', ypar = 5)$y / 5
  }
  
  points(seq(0, 6, 0.01), sapply(seq(0, 6, 0.01), muPq_perm_Halstead_uncertainty_DRC),
         pch = 5, col=4, cex = 0.5)

#LW1949 analysis
perm.lw1949 = dataprep(dose = perm.sum$perm.c, ntot = perm.sum$perm.total, nfx = perm.sum$perm.d)
perm.lw1949$cbitpfx <- constrain(perm.lw1949$bitpfx, probit(c(5e-04, 0.9995)))
perm.lwmod = lm(cbitpfx ~ log10dose, data = perm.lw1949[perm.lw1949$LWkeep, ])

  lines(seq(0, 6, 0.01), 
        pnorm(predict(perm.lwmod, newdata = data.frame(log10dose = log10(seq(0, 6, 0.01))), 
                      interval = 'confidence', level = 0.95)[,1]), lty = 2, col = 4)  
  
  lines(seq(0, 6, 0.01), 
        pnorm(predict(perm.lwmod, newdata = data.frame(log10dose = log10(seq(0, 6, 0.01))), 
                      interval = 'confidence', level = 0.95)[,2]), lty = 3, col = 4) 
  lines(seq(0, 6, 0.01), 
        pnorm(predict(perm.lwmod, newdata = data.frame(log10dose = log10(seq(0, 6, 0.01))), 
                      interval = 'confidence', level = 0.95)[,3]), lty = 3, col = 4) 

muPq_perm_Halstead_uncertainty = function(In){
  if(In==0) mup = 0 else{
    init = predict(perm.lwmod, newdata = data.frame(log10dose = log10(In)), se.fit = T)
    mup = pnorm(rnorm(1, init$fit, init$se.fit))# daily
  }
  return(mup)
}

points(seq(0, 6, 0.01), sapply(seq(0, 6, 0.01), muPq_perm_Halstead_uncertainty),
       pch = 5, col=6, cex = 0.5)

#Insecticide info from Halstead et al chemosphere paper ################
#calculate decay rates (k) from hydrolysis half lives for each chemical in Halstead 2015 (from table S1 and pmep.cce.cornell.edu) 
  mal.k = -log(0.5)/6.2         #from table s1 and in agreement of "less than 1 week in raw river water" from Cornell
  chlor.k = -log(0.5)/25.5      #from table S1; within the range of reported half life from cornell
  terb.k = -log(0.5)/6.5        #from table s1; in agreement with Cornell estimate of 5.5 days at pH of 7; degrades into formaldehyde
  lamcy.k = -log(0.5)/1         #very fast; 0 in table S1; "Not expected to be prevalent in surface waters" according to cornell website
  esfen.k = -log(0.5)/10        #cornell 4-15 days half life in water
  perm.k = -log(0.5)/2          #cornell "half life of less than 2.5 days"
  
  #median EECs from Halstead 2015
  med.mal = 0.778
  med.chlor = 5.810
  med.terb = 1.435
  med.lamcy = 0.649
  med.esfen = 0.311
  med.perm = 1.420
  
  #suggested application intervals for each agrochemical from Halstead 2015
  mal.days = c(1,6)             #two applications 5 days apart
  chlor.days = c(1,11,21)       #three applications 10 days apart
  terb.days = c(1)              #single application
  lamcy.days = seq(1, 61, by=4) #16 applications 4 days apart
  esfen.days = seq(1, 21, by=5) #5 applications 5 days apart
  perm.days = seq(1, 36, by=5)  #8 applications 5 days apart
  

#Character vector of objects to keep from this script ###############
  keep.hal15.muP = c('mal.sum', 'chlor.sum', 'terb.sum', 'esfen.sum', 'lamcy.sum', 'perm.sum',
                     'mal.mod', 'chlor.mod', 'terb.mod', 'esfen.mod', 'lamcy.mod', 'perm.mod',
                     'chlor.lwmod', 'terb.lwmod', 'esfen.lwmod', 'lamcy.lwmod', 'perm.lwmod',
                     'par.mal', 'par.chlor', 'par.terb', 'par.esfen', 'par.lamcy', 'par.perm',
                     'muPq_mal_Halstead_uncertainty', 'muPq_chlor_Halstead_uncertainty', 'muPq_terb_Halstead_uncertainty',
                     'muPq_esfen_Halstead_uncertainty', 'muPq_lamcy_Halstead_uncertainty', 'muPq_perm_Halstead_uncertainty')