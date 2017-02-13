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

ll4 = function(hi,lo,slp,lc,x){
  lo + ((hi-lo)/(1+exp(slp*(log(x)-lc))))
}

#Cercarial mortality (S. mansoni) from Tchounwou 1992 ##################
cerc = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/tchounwou92.csv')
  cerc.t = subset(cerc, chem == 'malathion')
  cerc.t$dead = round(cerc.t$dead)
  time = seq(0,25,0.1)

#Tchounwou Data plotted ############
plot(cerc.t$time_hrs[cerc.t$conc==0], cerc.t$surv[cerc.t$conc==0]/100, pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(cerc.t$conc[cerc.t$chem == 'malathion']))){
    points(cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]], 
           cerc.t$surv[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]]/100, pch=16,
           col = i)
  }
legend('topright', title = 'Mal(ppm)', legend = c(0,50,100,150,200,250), pch = c(17,rep(16,4)),
       col = c(1:6), cex=0.7)

#Fit models to survival curves and Add lines of best fit to plot ######################################
#Control ********************************************************************************************
  cerc.ctrl = subset(cerc.t, conc == 0)
  ctrl.mod = drm(alive/total ~ time_hrs, total, data = cerc.ctrl, 
                 type = 'binomial', fct = LL2.2())

  lines(time, ll4(1,0,ctrl.mod$coefficients[1], ctrl.mod$coefficients[2], time), lty=2)
  
#50 ppm ********************************************************************************************  
  cerc.mal50 = subset(cerc.t, conc == 50000)
  mal50.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal50, 
                  type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal50.mod$coefficients[1], mal50.mod$coefficients[2], time), lty=2, col=2)
  
#100ppm ********************************************************************************************
  cerc.mal100 = subset(cerc.t, conc == 100000)
  mal100.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal100, 
                   type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal100.mod$coefficients[1], mal100.mod$coefficients[2], time), lty=2, col=3)
  
#150 ppm ********************************************************************************************
  cerc.mal150 = subset(cerc.t, conc == 150000)
  mal150.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal150, 
                   type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal150.mod$coefficients[1], mal150.mod$coefficients[2], time), lty=2, col=4)
  
#200 ppm ********************************************************************************************
  cerc.mal200 = subset(cerc.t, conc == 200000)
  mal200.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal200, 
                   type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal200.mod$coefficients[1], mal200.mod$coefficients[2], time), lty=2, col=5)
  
#250ppm ********************************************************************************************
  cerc.mal250 = subset(cerc.t, conc == 250000)
  mal250.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal250, 
                   type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal250.mod$coefficients[1], mal250.mod$coefficients[2], time), lty=2, col=6)

#Data frame of LL.2 parameters across maltahion concentrations ##################    
parms.df = data.frame(mal = c(0,50,100,150,200,250),
                      e = c(coef(ctrl.mod)[2], coef(mal50.mod)[2], coef(mal100.mod)[2],
                            coef(mal150.mod)[2], coef(mal200.mod)[2], coef(mal250.mod)[2]),
                      e.se = c(summary(ctrl.mod)$coefficients[2,2], summary(mal50.mod)$coefficients[2,2],
                               summary(mal100.mod)$coefficients[2,2], summary(mal150.mod)$coefficients[2,2],
                               summary(mal200.mod)$coefficients[2,2], summary(mal250.mod)$coefficients[2,2]),
                      b = c(coef(ctrl.mod)[1], coef(mal50.mod)[1], coef(mal100.mod)[1],
                            coef(mal150.mod)[1], coef(mal200.mod)[1], coef(mal250.mod)[1]),
                      b.se = c(summary(ctrl.mod)$coefficients[1,2], summary(mal50.mod)$coefficients[1,2],
                               summary(mal100.mod)$coefficients[1,2], summary(mal150.mod)$coefficients[1,2],
                               summary(mal200.mod)$coefficients[1,2], summary(mal250.mod)$coefficients[1,2]))

  plot(parms.df$mal, parms.df$e, pch = 16, xlab = 'malathion (ppb)', ylab = 'LL.2 Parameters',
       ylim = c(-4, 8))
    points(parms.df$mal, parms.df$b, pch = 17, col=2)
    for(i in 1:length(parms.df$mal)){
      segments(x0 = parms.df$mal[i], y0 = parms.df$e[i] + parms.df$e.se[i],
               x1 = parms.df$mal[i], y1 = parms.df$e[i] - parms.df$e.se[i])
      segments(x0 = parms.df$mal[i], y0 = parms.df$b[i] + parms.df$b.se[i],
               x1 = parms.df$mal[i], y1 = parms.df$b[i] - parms.df$b.se[i], col=2)
    }
    
  e.mod = lm(e ~ mal, weights = e.se^-1, data = parms.df) 
  b.mod = lm(b ~ mal, weights = b.se^-1, data = parms.df) 
  b.mod2 = nls(b ~ a*exp(c*mal), start  = list(a = 2.6, c=0.002), weights = b.se^-1, data = parms.df)
    #Linear model has much lower AIC, and is somewhat more conservative so let's go with it

  mod.df= data.frame(mal = c(0:250),
                     pred.e = 0,
                     pred.e.se = 0,
                     pred.b = 0,
                     pred.b.se = 0)
  
  mod.df[,2:3] = predict(e.mod, newdata = mod.df, se.fit = TRUE)[1:2]
  mod.df[,4:5] = predict(b.mod, newdata = mod.df, se.fit = TRUE)[1:2]
    
  lines(mod.df$mal, mod.df$pred.e, lty = 2)
    lines(mod.df$mal, mod.df$pred.e + 1.96*mod.df$pred.e.se, lty = 3)
    lines(mod.df$mal, mod.df$pred.e - 1.96*mod.df$pred.e.se, lty = 3)
    
  lines(mod.df$mal, mod.df$pred.b, lty = 2, col=2)
    lines(mod.df$mal, mod.df$pred.b + 1.96*mod.df$pred.b.se, lty = 3, col=2)
    lines(mod.df$mal, mod.df$pred.b - 1.96*mod.df$pred.b.se, lty = 3, col=2)

#Function to estimate survival curve as function of q #####################################
pred.fx = function(In){
  e = as.numeric(predict(e.mod, newdata = data.frame(mal = In/1000), se.fit = TRUE)[1:2])
  b = as.numeric(predict(b.mod, newdata = data.frame(mal = In/1000), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t)-e.use))))}, 
                  lower=0, upper=24)[1]$value
  auc
}  

#Final function and keep vector
piC.tch92_mal_unc = function(In){
  piC = pred.fx(In) / pred.fx(0)
  if(piC > 1) piC = 1
  else(return(piC))
}  

keep.tch92.beq = c('piC.tch92_mal_unc', 'pred.fx', 'e.mod', 'b.mod')  