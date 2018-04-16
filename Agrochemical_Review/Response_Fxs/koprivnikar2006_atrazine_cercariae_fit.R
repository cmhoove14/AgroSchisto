#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(drc)

L.3.fx = function(t, lc50 = lc50, slp = slp){
    1 / (1+exp(slp*log(t / lc50)))
}  

#Cercarial mortality (E. trivolvis) from Koprivnikar 2006 ##################
kop.c = read.csv('Agrochemical_Review/Response_Fxs/Data/Koprivnikar2006.csv')

  kop.c[,5:8] = kop.c[,5:8]/100 #Convert survival measures to proportions
  time = c(0:25)
  kop.c$total = 5 #reported cercariae per treatment is 5-17, 5 per group represents conservative estimate
  
kop.cc = subset(kop.c, chem == 'control')  

kop.mod = drm(surv ~ time_hrs, conc, weights = surv.se^-1, data = kop.c, type = 'binomial', 
              fct = LL.2(names = c('b', 'e'), 
                         fixed = c(NA, NA)))

#Get estimate of cercariae-hrs for each concentration    
kop.atr.aucs.cer = as.numeric()

  kop.atr.aucs.cer[1] = integrate(f = L.3.fx, lc50 = kop.mod$coefficients[4], slp = kop.mod$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in c(2:3)){
    kop.atr.aucs.cer[j] = integrate(f = L.3.fx, lc50 = kop.mod$coefficients[j+3], slp = kop.mod$coefficients[j],
                                    lower=0, upper=24)[1]$value  
  }
    
#Create data frame with parameter values and atrazine concentrations #######################
kopc.df = data.frame(atrazine = c(0,20,200),
                     logatr = log(c(0,20,200)+1),
                     e = c(summary(kop.mod)$coefficients[c(4:6),1]),
                     e.se = c(summary(kop.mod)$coefficients[c(4:6),2]),
                     b = c(summary(kop.mod)$coefficients[c(1:3),1]),
                     b.se = c(summary(kop.mod)$coefficients[c(1:3),2]))

#parameters as function of atrazine  
kop.atr.e.lin = lm(e ~ atrazine, weights = e.se^-1, data = kopc.df)

  kop.atr.e.lin.pred = function(atr){
    predict(kop.atr.e.lin, newdata = data.frame(atrazine = atr), 
            interval = 'confidence', level = 0.95)
  }

kop.atr.e.exp = lm(e ~ logatr, weights = e.se^-1, data = kopc.df)

  kop.atr.e.exp.pred = function(atr){
    predict(kop.atr.e.exp, newdata = data.frame(logatr = log(atr+1)), 
            interval = 'confidence', level = 0.95)
  }

  #AIC(kop.atr.e.lin, kop.atr.e.exp) #linear fits better

#slope model  
kop.atr.b.lin = lm(b ~ atrazine, weights = b.se^-1, data = kopc.df) 

  kop.atr.b.lin.pred = function(atr){
    predict(kop.atr.b.lin, newdata = data.frame(atrazine = atr), 
            interval = 'confidence', level = 0.95)
  }

kop.atr.b.exp = lm(b ~ logatr, weights = b.se^-1, data = kopc.df) 

  kop.atr.b.exp.pred = function(atr){
    predict(kop.atr.b.exp, newdata = data.frame(logatr = log(atr+1)), 
            interval = 'confidence', level = 0.95)
  }
  
  #AIC(kop.atr.b.lin, kop.atr.b.exp) #exponential fits better
  
#Create function to generate d-r function with linear fit to lc50 parameter#####################
auc.kop.atr0 = function(He){
  e0 = as.numeric(predict(kop.atr.e.lin, newdata = data.frame(atrazine = He), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(kop.atr.b.exp, newdata = data.frame(logatr = log(He+1)), se.fit = TRUE)[1:2])
      
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  auc0
}
    
piC_kop_atr_unc = function(He){
  e = as.numeric(predict(kop.atr.e.lin, newdata = data.frame(atrazine = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(kop.atr.b.exp, newdata = data.frame(logatr = log(He+1)), se.fit = TRUE)[1:2])
      
  e.use = rnorm(1, e[1], e[2])
        while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
        while(b.use < 0) b.use = rnorm(1, e[1], e[2])
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, stop.on.error = FALSE)[1]$value
  piC = auc/auc.kop.atr0(0) 

  return(piC)
}  
    
keep.kop06.beq = c('auc.kop.atr0', 'piC_kop_atr_unc', 'kop.atr.e.lin', 'kop.atr.b.exp')        