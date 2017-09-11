#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#This code adapted from Arathi's "schisto_halstead_2pops_mda_bouncebackRate_cloudParam.R"

##################################################################################################

require(deSolve)
require(graphics)
require(ggplot2)

source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")

outputfile<-"C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/shortlist_51.csv"
shortlist_first100<-as.data.frame(read.csv(file=outputfile, header=TRUE, sep=","))
sort_ind<-order(shortlist_first100$negLL, decreasing=FALSE)   #Sorted by max of neg LL
shortlist_first100<-shortlist_first100[sort_ind,]
  shortlist_first100 = shortlist_first100[c(1:100),] #Get rid of other parameter estimates
  
params = parameters_2pops_mda_Chris1 
params['beta'] = shortlist_first100$beta[1]
params['lamda'] = shortlist_first100$lamda.twa[1]
 
#Mating probability function
fx<-function(x, mean.worm, clump){
  alpha<-(mean.worm)/(clump+mean.worm)
  (1-cos(x)) / ( (1+(alpha*cos(x)))^(1+clump) )
}

phi_Wk<-function(W, k){
  alpha<-W/(W+k)
  1-( (1-alpha)^(k+1) * (integrate(fx, 0, 2*pi, W, k)$value) /(2*pi)  )
} 
 
#Parameter values 
    f_N<-params["f_N"]
    phi_N<-params["phi_N"]
    z<-params["z"]
    mu_N<-params["mu_N"]
    sigma<-params["sigma"]
    mu_I<-params["mu_I"]
    mu_W<-params["mu_W"]
    m<-params["m"]
    H<-params["H"]
    mu_H<-params["mu_H"]
    muPq<-params["muPq"]
    beta<-params["beta"] 
    lamda<-params["lamda"]
    k<-params['k']
    
#Reff function
  getReff<-function(W, k = 0.08){ 
    
    k1 = sigma/(mu_N+mu_I)
    k2 = beta*0.5*W*H*phi_Wk(W,k)/(mu_N+sigma)
    
    Num_1 = lamda*k1*k2*phi_N
    Num_2 = ( (1+k2)*f_N ) - ( mu_N + (beta*0.5*H*W*phi_Wk(W,k) ) )
    Den = ( 1+k2+(k1*k2) ) * ( 1+k2 ) * ( mu_W+mu_H )* W* f_N
    
    Reff = as.numeric(Num_1*Num_2/Den)
    
    k2_nophi = beta*0.5*W*H/(mu_N+sigma)
    
    Num_1_nophi = lamda*k1*k2_nophi*phi_N
    Num_2_nophi = ( (1+k2_nophi)*f_N ) - ( mu_N + (beta*0.5*H*W ) )
    Den_nophi = ( 1+k2_nophi+(k1*k2_nophi) ) * ( 1+k2_nophi ) * ( mu_W+mu_H )* W* f_N
    
    Reff_noMatingFn = as.numeric(Num_1_nophi*Num_2_nophi/Den_nophi)
    
    return(c(Reff, Reff_noMatingFn))
  }
  
#Vectors to fill
  W.seq = c(seq(from=0.01, to=1, by = 0.005), seq(from=2, to=200, by = 0.05))
  Reff.seq = as.numeric()
  r0.seq = as.numeric()
  dwdt.seq = as.numeric()
  dwdt.seq2 = as.numeric()

#Fill vectors across range of W values    
for(w in 1:length(W.seq)){
    W = W.seq[w]
    Reff.seq[w] = getReff(W, k)[1]
    r0.seq[w] = getReff(W, k)[2]

    dwdt.seq[w] = (Reff.seq[w]-1)*W*(mu_W+mu_H) 
    dwdt.seq2[w] = (r0.seq[w]-1)*W*(mu_W+mu_H) 
    
}
  
opar<-par()
    
#plot Reff profile with log(W) ##########
windows(width = 19, height = 11)    
plot(log(W.seq), Reff.seq, type = 'l', lwd = 3, 
     ylim = c(0,max(Reff.seq)), xlim = log(c(min(W.seq)+0.003, max(W.seq))),
     xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')  
  abline(h = 1, lty=3, lwd = 2)
  #abline(h = max(r0.seq), lty = 3, col = 2)
  lines(log(W.seq), r0.seq, lty = 2, lwd = 3)
  
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), #col = 3,
             y0 = 1, y1 = -1, lty = 2, lwd = 2)
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4]), #col = 3,
             y0 = 1, y1 = -1, lty = 2, lwd = 2)
    segments(x0 = log(W.seq[which(Reff.seq == max(Reff.seq))]),
             x1 = log(W.seq[which(Reff.seq == max(Reff.seq))]), #col = 3,
             y0 = max(Reff.seq), y1 = -1, lty = 2, lwd = 2)
    
    axis(1, at = c(log(W.seq[which(round(Reff.seq, digits = 2) <= 1.0)][1]), 
                   log(W.seq[which(Reff.seq == max(Reff.seq))]),
                   log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4])),
         labels = c(expression(paste(italic('W'[bp]), sep = '')),
                    expression(paste(italic(' W'[peak]), sep = '')),
                    expression(paste(italic('    W'[eq]), sep = ''))), #col = 3,
         lty = 2, cex.axis = 1.2)
    axis(1, at = c(-4,0,4), mgp = c(3,1.2,0), cex.axis = 1.2)
    axis(2, at = max(r0.seq), labels = expression(paste(italic('R'[0]), sep = '')), #col = 2, 
         lty=2, las = 2, mgp = c(3,0.75,0), cex.axis = 1.2)
    axis(2, at = c(0,1,2,3), labels = c('0', '1', '2', '3'), cex.axis = 1.2)
    
    mtext(side = 2, text = expression(paste('R'[eff], sep = '')), line = 2.5, cex = 1.4)
    mtext(side = 1, text = expression(paste('     log (', italic('W'), ')', sep = '')),
          line = 3, cex = 1.4)
    
#plot Reff profile with W ##########
  windows(width = 19, height = 11)    
    plot(W.seq, Reff.seq, type = 'l', lwd = 3, 
         ylim = c(0,max(Reff.seq)), xlim = c(0, max(W.seq)),
         xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')  
    abline(h = 1, lty=3, lwd = 2)
    #abline(h = max(r0.seq), lty = 3, col = 2)
    lines(W.seq, r0.seq, lty = 2, lwd = 3)
    
    segments(x0 = W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1],
             x1 = W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1], #col = 3,
             y0 = 1, y1 = -1, lty = 2, lwd = 2)
    segments(x0 = W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4],
             x1 = W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4], #col = 3,
             y0 = 1, y1 = -1, lty = 2, lwd = 2)
    segments(x0 = W.seq[which(Reff.seq == max(Reff.seq))],
             x1 = W.seq[which(Reff.seq == max(Reff.seq))], #col = 3,
             y0 = max(Reff.seq), y1 = -1, lty = 2, lwd = 2)
    
    axis(1, at = c(W.seq[which(round(Reff.seq, digits = 1) <= 1.0)][1], 
                   W.seq[which(Reff.seq == max(Reff.seq))],
                   W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4]),
         labels = c(expression(paste(italic('W'[bp]), sep = '')),
                    expression(paste(italic(' W'[peak]), sep = '')),
                    expression(paste(italic('    W'[eq]), sep = ''))), #col = 3,
         lty = 2, cex.axis = 1.2)
    axis(1, at = c(0,20,40,60,80,100), mgp = c(3,1.2,0), cex.axis = 1.2)
    axis(2, at = max(r0.seq), labels = expression(paste(italic('R'[0]), sep = '')), #col = 2, 
         lty=2, las = 2, mgp = c(3,0.75,0), cex.axis = 1.2)
    axis(2, at = c(0:8), labels = c('0', '1', '2', '3','4','5','6','7','8'), cex.axis = 1.2)
    
    mtext(side = 2, text = expression(paste('R'[eff], sep = '')), line = 2.5, cex = 1.4)
    mtext(side = 1, text = expression(paste(italic('W'), sep = '')),
          line = 3, cex = 1.4)
    
    
#Increase snail mortality rate and observe effect on Reff profile ##############
  mu_N<-params["mu_N"]*2
  
#Vectors to fill
  W.seq = c(seq(from=0.01, to=1, by = 0.005), seq(from=2, to=200, by = 0.05))
  Reff.seq.mun2 = as.numeric()
  r0.seq.mun2 = as.numeric()
  dwdt.seq.mun2 = as.numeric()
  dwdt.seq2.mun2 = as.numeric()
  
#Fill vectors across range of W values    
  for(w in 1:length(W.seq)){
    W = W.seq[w]
    Reff.seq.mun2[w] = getReff(W, k)[1]
    r0.seq.mun2[w] = getReff(W, k)[2]
    
    dwdt.seq.mun2[w] = (Reff.seq.mun2[w]-1)*W*(mu_W+mu_H) 
    dwdt.seq2.mun2[w] = (r0.seq.mun2[w]-1)*W*(mu_W+mu_H) 
    
  }
  
  lines(W.seq, Reff.seq.mun2, col = 2, lwd = 2)
  
#Increase adult worm mortality rate and observe effect on Reff profile ##############
  mu_N<-params["mu_N"]   #Reset snail mortality to baseline
  mu_W<-params["mu_W"]*2 #double adult worm mortality
  
  #Vectors to fill
  W.seq = c(seq(from=0.01, to=1, by = 0.005), seq(from=2, to=200, by = 0.05))
  Reff.seq.muw2 = as.numeric()
  r0.seq.muw2 = as.numeric()
  dwdt.seq.muw2 = as.numeric()
  dwdt.seq2.muw2 = as.numeric()
  
  #Fill vectors across range of W values    
  for(w in 1:length(W.seq)){
    W = W.seq[w]
    Reff.seq.muw2[w] = getReff(W, k)[1]
    r0.seq.muw2[w] = getReff(W, k)[2]
    
    dwdt.seq.muw2[w] = (Reff.seq.muw2[w]-1)*W*(mu_W+mu_H) 
    dwdt.seq2.muw2[w] = (r0.seq.muw2[w]-1)*W*(mu_W+mu_H) 
    
  }
  
  lines(W.seq, Reff.seq.muw2, col = 3, lwd = 2)
    
#Save plot as high res tiff ##################
tiff("Elimination_Feasibility/plots/PLoS_Figs/Fig1.tiff", height = 11, width = 19.05, units = 'cm', 
     compression = "lzw", res = 300)
plot(log(W.seq), Reff.seq, type = 'l', lwd = 3, 
     ylim = c(0,2), xlim = log(c(min(W.seq)+0.003, max(W.seq))),
     xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')  
abline(h = 1, lty=3, lwd = 2)
#abline(h = max(r0.seq), lty = 3, col = 2)
lines(log(W.seq), r0.seq, lty = 2, lwd = 3)

segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
         x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), #col = 3,
         y0 = 1, y1 = -1, lty = 2, lwd = 2)
segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4]),
         x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4]), #col = 3,
         y0 = 1, y1 = -1, lty = 2, lwd = 2)
segments(x0 = log(W.seq[which(Reff.seq == max(Reff.seq))]),
         x1 = log(W.seq[which(Reff.seq == max(Reff.seq))]), #col = 3,
         y0 = max(Reff.seq), y1 = -1, lty = 2, lwd = 2)

axis(1, at = c(log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), 
               log(W.seq[which(Reff.seq == max(Reff.seq))]),
               log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4])),
     labels = c(expression(paste(italic('W'[bp]), sep = '')),
                expression(paste(italic(' W'[peak]), sep = '')),
                expression(paste(italic('    W'[eq]), sep = ''))), #col = 3,
     lty = 2,  mgp = c(3,1,0), cex.axis = 1.2)
axis(1, at = c(-4,0,4), mgp = c(3,1.2,0), cex.axis = 1.2)
axis(2, at = max(r0.seq), labels = expression(paste(italic('R'[0]), sep = '')), #col = 2, 
     lty=2, las = 2, mgp = c(3,0.75,0), cex.axis = 1.2)
axis(2, at = c(0,1,2,3), labels = c('0', '1', '2', '3'), cex.axis = 1.2)

mtext(side = 2, text = expression(paste('R'[eff], sep = '')), line = 2.5, cex = 1.4)
mtext(side = 1, text = expression(paste('     log (', italic('W'), ')', sep = '')),
      line = 3, cex = 1.4)
dev.off()

