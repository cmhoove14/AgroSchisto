#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

atr = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/rohr08_atr.csv')

plot(atr$log_conc, atr$surv, xlab = 'log+1 Atrazine (ppb)', ylab = 'Prop alive @14-18 hrs', 
     xlim = c(0,8), ylim = c(0.25,0.65), pch = 16)
  for(i in 1:length(atr$log_conc)){
    segments(x0 = atr$log_conc[i], y0 = 1-atr$hi[i],
             x1 = atr$log_conc[i], y1 = 1-atr$lo[i])
  }

atr.lin = lm(surv ~ log_conc, weights = st_err^-1, data = atr)
rohr.atr.fx = function(He){
  heu = log(He+1)
  predict(atr.lin, newdata = data.frame(log_conc = heu), interval = 'confidence', level = 0.95)
}

  lines(log(seq(0,2000,20)+1), sapply(seq(0,2000,20), rohr.atr.fx)[1,], lty = 2)
  lines(log(seq(0,2000,20)+1), sapply(seq(0,2000,20), rohr.atr.fx)[2,], lty = 3)
  lines(log(seq(0,2000,20)+1), sapply(seq(0,2000,20), rohr.atr.fx)[3,], lty = 3)
      
#Final function estimates relative mortality to He=0 and is interpreted as pi_C parameter      
  piC.atr.rohr08.lin = function(He){
      ts = predict(atr.lin, newdata = data.frame(log_conc = log(He+1)), se.fit = TRUE)[1:2]
      piC = rnorm(1, ts$fit, ts$se.fit) / atr.lin$coefficients[1]
    
    return(piC)
  } 
  
  
#Plot to see how it does
  plot(atr$conc, atr$surv / atr$surv[1], pch = 16, xlim = c(0,2000), ylim = c(0,1),
       ylab = 'rel prop alive @ 14-18 hrs', xlab = 'Atrazine (ppb)',
       main = 'Sample output of cercarial mortality')
    points(seq(0, 2000, 5), sapply(seq(0, 2000, 5), piC.atr.rohr08.lin), pch = 17, col=3, cex = 0.5) 
    legend('bottomleft', legend = c('Observed', 'Modeled'), 
           pch = c(16,17), col = c(1,3), cex = 0.7, bty = 'n')

#Vector of items to keep
  keep.atr.rohr08 = c('atr', 'atr.lin', 'piC.atr.rohr08.lin')