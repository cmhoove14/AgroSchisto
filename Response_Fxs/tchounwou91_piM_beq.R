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
  lo + ((hi-lo)/(1+exp(slp*(log(x/lc)))))
}
L.3.fx = function(t, lc50 = lc50, slp = slp){
  1 / (1+exp(slp*log(t / lc50)))
}
  time = seq(0,25,0.1)

#miracidial mortality (S. mansoni) from Tchounwou 1991 ####################################
mir = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Miracidial Mortality/tchounwou08_miracidia.csv')
  mir.mal = subset(mir, chem == 'mal')

#Tchounwou Data plotted ############
tch91.piM.mal<-drm(alive/total ~ time_hrs, conc, weights = total, data = mir.mal, type = 'binomial', 
                   fct = LL.4(names = c('b', 'c', 'd', 'e'),
                              fixed = c(NA, 0, 1, NA)))
  summary(tch91.piM.mal)
  plot(tch91.piM.mal) 

plot(mir.mal$time_hrs[mir.mal$conc==0], mir.mal$alive[mir.mal$conc==0]/mir.mal$total[mir.mal$conc==0], 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tch91.piM.mal, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(mir.mal$conc))){
    points(mir.mal$time_hrs[mir.mal$conc==unique(mir.mal$conc)[i]], 
           mir.mal$alive[mir.mal$conc==unique(mir.mal$conc)[i]] /
           mir.mal$total[mir.mal$conc==unique(mir.mal$conc)[i]], pch=16,
           col = i)
    lines(time, predict(tch91.piM.mal, 
                        data.frame(time_hrs=time, conc = unique(mir.mal$conc)[i])),
          lty = 2, col = i)
  }

title('malathion toxicity to cercariae')
legend('topright', legend = c('control', unique(mir.mal$conc)[-1]), 
       pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')  

#Get estimate of miracidia-hrs for each concentration    
tch91.mal.pim.aucs = as.numeric()
  
  for(j in 1:length(unique(mir.mal$conc))){
    fx = function(t){
      predict(tch91.piM.mal, newdata = data.frame(time_hrs = t, conc = unique(mir.mal$conc)[j]))
    }
    tch91.mal.pim.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }


#Create data frame with parameter values and malathion concentrations #######################
mirp.df = data.frame(mal = c(0,30,60,90,120,150)*1000,
                     logmal = log(c(0,30,60,90,120,150)+1*1000),
                     e = summary(tch91.piM.mal)$coefficients[c(7:12), 1],
                     e.se = summary(tch91.piM.mal)$coefficients[c(7:12), 2],
                     b = summary(tch91.piM.mal)$coefficients[c(1:6), 1],
                     b.se = summary(tch91.piM.mal)$coefficients[c(1:6), 2])

plot(mirp.df$mal, mirp.df$e, pch = 16, xlab = 'malathion (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0, 7))
  points(mirp.df$mal, mirp.df$b, pch = 17, col=2)
    for(i in 1:length(mirp.df$mal)){
      segments(x0 = mirp.df$mal[i], y0 = mirp.df$e[i] + mirp.df$e.se[i],
               x1 = mirp.df$mal[i], y1 = mirp.df$e[i] - mirp.df$e.se[i])
      segments(x0 = mirp.df$mal[i], y0 = mirp.df$b[i] + mirp.df$b.se[i],
               x1 = mirp.df$mal[i], y1 = mirp.df$b[i] - mirp.df$b.se[i], col=2)
    }

tch91.em.mod = lm(e ~ mal, weights = e.se^-1, data = mirp.df) 
  tch91.em.fx = function(mal){
    predict(tch91.em.mod, newdata = data.frame(mal = mal), interval = 'confidence', level = 0.95)
  }
    lines(seq(0,max(mirp.df$mal)+1000,100), sapply(seq(0,max(mirp.df$mal)+1000,100), tch91.em.fx)[1,],
          lty = 2)
    lines(seq(0,max(mirp.df$mal)+1000,100), sapply(seq(0,max(mirp.df$mal)+1000,100), tch91.em.fx)[2,],
          lty = 3)
    lines(seq(0,max(mirp.df$mal)+1000,100), sapply(seq(0,max(mirp.df$mal)+1000,100), tch91.em.fx)[3,],
          lty = 3)
#not trying exponential because linear is clearly the best fit
    
tch91.bm.mod = lm(b ~ mal, weights = b.se^-1, data = mirp.df) 
  tch91.bm.fx = function(mal){
    predict(tch91.bm.mod, newdata = data.frame(mal = mal), interval = 'confidence', level = 0.95)
  }
    lines(seq(0,max(mirp.df$mal)+1000,100), sapply(seq(0,max(mirp.df$mal)+1000,100), tch91.bm.fx)[1,],
          lty = 2, col = 2)
    lines(seq(0,max(mirp.df$mal)+1000,100), sapply(seq(0,max(mirp.df$mal)+1000,100), tch91.bm.fx)[2,],
          lty = 3, col = 2)
    lines(seq(0,max(mirp.df$mal)+1000,100), sapply(seq(0,max(mirp.df$mal)+1000,100), tch91.bm.fx)[3,],
          lty = 3, col = 2)

mir.auc.mal = data.frame(mal = mirp.df$mal,
                         piM = c(tch91.mal.pim.aucs)/tch91.mal.pim.aucs[1])    
#Create function to generate d-r function #####################
piM.tch91_mal_unc = function(In){
  if(In == 0) piM = 1 else{
  e = as.numeric(predict(tch91.em.mod, newdata = data.frame(mal = In), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.bm.mod, newdata = data.frame(mal = In), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  while(b.use < 0) b.use = rnorm(1, e[1], e[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
  piM = auc/tch91.mal.pim.aucs[1]
  if(piM > 1) piM = 1
  }
  return(piM)
}  
  
keep.tch91.beq = c('tch91.mal.pim.aucs', 'piM.tch91_mal_unc', 'L.3.fx', 'mirp.df',
                   'tch91.em.mod', 'tch91.bm.mod', 'mir.auc.mal')

#Qualitative model validation ###############
#Regenerate plot of observed data
plot(mir.mal$time_hrs[mir.mal$conc==0], mir.mal$surv[mir.mal$conc==0], pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(mir.mal$conc))){
    points(mir.mal$time_hrs[mir.mal$conc==unique(mir.mal$conc)[i]], 
           mir.mal$surv[mir.mal$conc==unique(mir.mal$conc)[i]], pch=16,
           col = i)
  }
  legend('topright', legend = c(0, 30, 60, 90, 120, 150), title = 'Malathion (ppm)',
         pch = c(17,rep(16,5)), col = c(1,2:6), cex = 0.8, bty = 'n')

#function to plot model predictions
predm.fx.plot = function(In, clr){
  e = as.numeric(predict(tch91.em.mod, newdata = data.frame(mal = In), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.bm.mod, newdata = data.frame(mal = In), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  while(b.use < 0) b.use = rnorm(1, e[1], e[2])
  
  lines(time, ll4(1,0, b.use, e.use, time), lty=2, col = clr)
}   

#plot model predictions
for(i in c(0,30,60,90,120,150)*1000){
  c = i/30000 + 1
  print(c)
  replicate(10, predm.fx.plot(In = i, clr = c))
}

#plot model output compared to observed points
plot(mirp.df$mal, tch91.mal.pim.aucs/tch91.mal.pim.aucs[1], pch = 16, ylim = c(0,1),
     xlab = 'Malathion (ppb)', ylab = expression(paste(pi[M])),
     main = 'Sample Output of miracidial mortality function')
  points(seq(0,max(mirp.df$mal)+1000,250), sapply(seq(0,max(mirp.df$mal)+1000,250), piM.tch91_mal_unc), 
         pch = 5, col=4, cex = 0.5)