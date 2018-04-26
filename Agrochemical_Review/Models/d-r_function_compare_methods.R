#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Due to difficulty reproducing d-r function from Abdel-Ghaffar 2016 paper, sought to establish another way to reproduce dose-response function from reported Lithfield and Wilcoxon 1949 parameters. Abdel-Ghaffar's slope and LC50 values did not reproduce a d-r function that matched their other reported data (LC0, LC10, LC90). Here using similar data from Bakry et al 2011, compared the reproduced litchfield and wilcoxon method with estimates from linear function fit to reported LC values. Results show that the mean mortality estimates are very similar at low concentrations that are within the expected environmental concentration ofeach chemical. In addition, the variance of the linear method is higher than that of the reproduced litchfield and wilcozxon method, implying that this method more conservative in addition to approximating the mean at low concentrations.

source("Agrochemical_Review/Models/litchfield_wilcoxon_get_b1_from_slope.R")
#Snail (H. duryi) toxicity ##########
mun.mal = data.frame(conc = c(0, .176, .480, .830, 1.760, 3.360)*1000,
                     mort = c(0, 0,   .10, .25, .50 , .90),
                     surv = 0)
  mun.mal$surv = 1 - mun.mal$mort
  mun.mal$ppm = mun.mal$conc/1000
  mun.mal$ppmlog10 = log10(mun.mal$ppm)
  mun.mal$probit = qnorm(mun.mal$mort)
  
  mun.mal.sub = subset(mun.mal, is.finite(probit))

#Create function based on LC values
  plot(mun.mal.sub$ppmlog10, mun.mal.sub$probit, pch = 16, xlab = "log10 concentration", ylab = "probit mortality (LC)")
    abline(lm(probit~ppmlog10, data = mun.mal.sub), lty = 2, col = 3)
    
  bak_mal_lin <- lm(probit~ppmlog10, data = mun.mal.sub)
    b.lin <- coef(bak_mal_lin)[1]
      b.lin.up <- confint(bak_mal_lin)[1,1]
      b.lin.lo <- confint(bak_mal_lin)[1,2]
    m.lin <- coef(bak_mal_lin)[2]
  
muNq_mal_Bakry11_lin <- function(In, b = b.lin){
  mun = pnorm(b + log10(In/1000)*m.lin)
  return(mun)
}  

lin_pred <- predict(bak_mal_lin, newdata = data.frame(ppmlog10 = log10(c(0:5000)/1000)), interval = "confidence", level = 0.95)

lin_pred_trnsfrm <- pnorm(lin_pred)

#Malathion reported LC50 and slope data #####################
  lc50.bak.mal.report = 1.760
  slp.bak.mal.report = 2.68
  b1.bak.mal = get_b1(slp.bak.mal.report)
  #get standard error from reported 95% CIs of lc50
      se.lc50.bak.mal = mean(c(log10(3.12/lc50.bak.mal.report), log10(lc50.bak.mal.report/0.99))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
muNq_mal_Bakry11_uncertainty = function(In){
  Ins = (In/1000) #Parameters based on ppm, data input as ppb
  lc50 = 10^(rnorm(1, log10(lc50.bak.mal.report), se.lc50.bak.mal)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm(b1.bak.mal * log10(Ins/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }
  
muNq_mal_Bakry11 = function(In, lc50){
  Ins = (In/1000) #Parameters based on ppm, data input as ppb

    mun = pnorm(b1.bak.mal * log10(Ins/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }

plot(mun.mal$conc, mun.mal$mort, pch = 16, 
     xlab = 'Malathion (ppb)', ylab = expression(italic(mu[N])), ylim = c(0,1), xlim = c(0,5000),
     main = 'Bakry et al 2011 malathion snail mortality')
  segments(x0 = 3120, y0 = 0.5, x1 = 990, y1 = 0.5)
  
  set.seed(43093)
  lines(c(0:5000), sapply(c(0:5000), muNq_mal_Bakry11, lc50 = lc50.bak.mal.report, simplify = TRUE), lty = 2, col = 2)
  lines(c(0:5000), sapply(c(0:5000), muNq_mal_Bakry11, lc50 = 3.12, simplify = TRUE), lty = 3, col = 2)
  lines(c(0:5000), sapply(c(0:5000), muNq_mal_Bakry11, lc50 = 0.99, simplify = TRUE), lty = 3, col = 2)
  
  lines(c(0:5000), sapply(c(0:5000), muNq_mal_Bakry11_lin, simplify = TRUE), lty = 2, col = 3)
  lines(c(0:5000), sapply(c(0:5000), muNq_mal_Bakry11_lin, b = b.lin.up, simplify = TRUE), lty = 3, col = 3)
  lines(c(0:5000), sapply(c(0:5000), muNq_mal_Bakry11_lin, b = b.lin.lo, simplify = TRUE), lty = 3, col = 3)
  
  lines(c(0:5000), lin_pred_trnsfrm[,1], lty = 2, col = 4)
  lines(c(0:5000), lin_pred_trnsfrm[,2], lty = 3, col = 4)
  lines(c(0:5000), lin_pred_trnsfrm[,3], lty = 3, col = 4)

  legend("bottomright", bty="n", lty = 2, col = c(2,3,4), legend = c("LW1949", "LC_ fit, intercept uncertainty", "LC_ fit, slope and intercept uncertainty"), cex = 0.75)

#keep vector
  keep.bak.mal = c('muNq_mal_Bakry11_uncertainty', 'lc50.bak.mal.report', 'se.lc50.bak.mal', 'b1.bak.mal')    
  
#Deltamethrin reported LC50 and slope data #############  
mun.del = data.frame(conc = c(0, .482, 1.21, 2.034, 4.82, 7.26)*1000,
                     mort = c(0, 0,   .10, .25, .50 , .90),
                     surv = 0)
  mun.del$surv = 1 - mun.del$mort
  mun.del$ppm = mun.del$conc/1000
  mun.del$ppmlog10 = log10(mun.del$ppm)
  mun.del$probit = qnorm(mun.del$mort)

  mun.del.sub = subset(mun.del, is.finite(probit))

#Create function based on LC values
plot(mun.del.sub$ppmlog10, mun.del.sub$probit, pch = 16, xlab = "log10 concentration", ylab = "probit mortality (LC)")
    abline(lm(probit~ppmlog10, data = mun.del.sub), lty = 2, col = 3)
    
  bak_del_lin <- lm(probit~ppmlog10, data = mun.del.sub)
    b.lin.del <- coef(bak_del_lin)[1]
      b.lin.up.del <- confint(bak_del_lin)[1,1]
      b.lin.lo.del <- confint(bak_del_lin)[1,2]
    m.lin.del <- coef(bak_del_lin)[2]
  
  muNq_del_Bakry11_lin <- function(In, b = b.lin.del){
    mun = pnorm(b + log10(In/1000)*m.lin.del)
    return(mun)
  }  

  lin_del <- predict(bak_del_lin, newdata = data.frame(ppmlog10 = log10(c(0:10000)/1000)), interval = "confidence", level = 0.95)
  
  lin_del_trnsfrm <- pnorm(lin_del)

#Function based on back transform from Litchfield and Wilcoxon  
  lc50.bak.del.report = 4.82
  slp.bak.del.report = 2.74
  b1.bak.del = get_b1(slp.bak.del.report)
  #get standard error from reported 95% CIs of lc50
    se.lc50.bak.del = mean(c(log10(7.7/lc50.bak.del.report), log10(lc50.bak.del.report/3.1))) / 1.96
  
  muNq_del_Bakry11_uncertainty = function(In){
    Ins = (In/1000) #Parameters based on ppm, data input as ppb
    lc50 = 10^(rnorm(1, log10(lc50.bak.del.report), se.lc50.bak.del)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm(b1.bak.del * log10(Ins/lc50))
    
    return(mun)
  }
  
muNq_del_Bakry11 = function(In, lc50){
  Ins = (In/1000) #Parameters based on ppm, data input as ppb

    mun = pnorm(b1.bak.del * log10(Ins/lc50)) #Estimate daily mortality (percent)

    return(mun)
}

plot(mun.del$conc, mun.del$mort, pch = 16, ylim = c(0,1), xlim = c(0,10000),
     xlab = 'Deltamethrin (ppb)', ylab = expression(italic(mu[N])), 
     main = 'Bakry et al 2011 deltamethrin snail mortality')
  segments(x0 = 3100, y0 = 0.5, x1 = 7700, y1 = 0.5)
  
  set.seed(43093)
  lines(c(0:10000), sapply(c(0:10000), muNq_del_Bakry11, lc50 = lc50.bak.del.report, simplify = TRUE), lty = 2, col = 2)
  lines(c(0:10000), sapply(c(0:10000), muNq_del_Bakry11, lc50 = 7.7, simplify = TRUE), lty = 3, col = 2)
  lines(c(0:10000), sapply(c(0:10000), muNq_del_Bakry11, lc50 = 3.1, simplify = TRUE), lty = 3, col = 2)
  
  lines(c(0:10000), sapply(c(0:10000), muNq_del_Bakry11_lin, simplify = TRUE), lty = 2, col = 3)
  lines(c(0:10000), sapply(c(0:10000), muNq_del_Bakry11_lin, b = b.lin.up.del, simplify = TRUE), lty = 3, col = 3)
  lines(c(0:10000), sapply(c(0:10000), muNq_del_Bakry11_lin, b = b.lin.lo.del, simplify = TRUE), lty = 3, col = 3)
  
  lines(c(0:10000), lin_del_trnsfrm[,1], lty = 2, col = 4)
  lines(c(0:10000), lin_del_trnsfrm[,2], lty = 3, col = 4)
  lines(c(0:10000), lin_del_trnsfrm[,3], lty = 3, col = 4)

  legend("bottomright", bty="n", lty = 2, col = c(2,3,4), legend = c("LW1949", "LC_ fit, intercept uncertainty", "LC_ fit, slope and intercept uncertainty"), cex = 0.75)
