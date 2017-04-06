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
sort_ind<-order(shortlist_first100$R0, decreasing=FALSE)
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
  
for(w in 1:length(W.seq)){
    W = W.seq[w]
    Reff.seq[w] = getReff(W, k)[1]
    r0.seq[w] = getReff(W, k)[2]

    dwdt.seq[w] = (Reff.seq[w]-1)*W*(mu_W+mu_H) 
    dwdt.seq2[w] = (r0.seq[w]-1)*W*(mu_W+mu_H) 
    
}
  
opar<-par()
    
#plot Reff profile ##########
windows(width = 8, height = 6)    
plot(log(W.seq), Reff.seq, type = 'l', lwd = 2, 
     ylim = c(0,2), xlim = log(c(min(W.seq)+0.003, max(W.seq))),
     xaxt = 'n', yaxt = 'n', ylab = '',
     xlab = expression(paste('log (', italic('W'), ')', sep = '')))  
  abline(h = 1, lty=3)
  #abline(h = max(r0.seq), lty = 3, col = 2)
  lines(log(W.seq), r0.seq, lty = 2, lwd = 2)
  
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), #col = 3,
             y0 = 1, y1 = -1, lty = 2)
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4]), #col = 3,
             y0 = 1, y1 = -1, lty = 2)
    segments(x0 = log(W.seq[which(Reff.seq == max(Reff.seq))]),
             x1 = log(W.seq[which(Reff.seq == max(Reff.seq))]), #col = 3,
             y0 = max(Reff.seq), y1 = -1, lty = 2)
    
    axis(1, at = c(log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), 
                   log(W.seq[which(Reff.seq == max(Reff.seq))]),
                   log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4])),
         labels = c(expression(paste(italic('W'[bp]), sep = '')),
                    expression(paste(italic(' W'[peak]), sep = '')),
                    expression(paste(italic('    W'[eq]), sep = ''))), #col = 3,
         lty = 2,  mgp = c(3,0.5,0), cex.axis = 0.8)
    axis(1, at = c(-4,0,4))
    axis(2, at = max(r0.seq), labels = expression(paste(italic('R'[0]), sep = '')), #col = 2, 
         lty=2, las = 2, mgp = c(3,0.75,0), cex.axis = 0.8)
    axis(2, at = c(0,1,2,3), labels = c('0', '1', '2', '3'))
    
    mtext(side = 2, text = expression(paste('R'[eff], sep = '')), line = 2)

#Plot BBR profile #######
    
windows(height = 8, width = 10)
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,
                2,2,2,2,2,2,2,2,2,2,
                2,3,3,3,3,2,2,2,2,2,
                2,3,3,3,3,2,2,2,2,2,
                2,3,3,3,3,2,2,2,2,2,
                2,2,2,2,2,2,2,2,2,2,
                2,2,2,2,2,2,2,2,2,2), 
                10, 10, byrow = TRUE))
  par(mar=c(4.1,4.7,1.6,1.1), cex.lab = 1.5)    
plot(log(W.seq), (dwdt.seq*365)/W.seq, type='l', lwd=2, ylim=c(-0.3, 0.3), xlim =c(-4.3, log(max(W.seq))),
     ylab="", xaxt = 'n', yaxt ='n', xlab = '')
    abline(h = 0, lty=3)
    lines(log(W.seq), (dwdt.seq2*365)/W.seq, lwd = 2, lty = 2)
    
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), #col = 3,
             y0 = 0, y1 = -1, lty = 2)
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4]), #col = 3,
             y0 = 0, y1 = -1, lty = 2)
    segments(x0 = log(W.seq[which(Reff.seq == max(Reff.seq))]),
             x1 = log(W.seq[which(Reff.seq == max(Reff.seq))]), #col = 3,
             y0 = max((dwdt.seq*365)/W.seq), y1 = -1, lty = 2)
    
    axis(1, at = c(log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), 
                   log(W.seq[which(Reff.seq == max(Reff.seq))]),
                   log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][4])),
         labels = c(expression(paste(italic('W'[bp]), sep = '')),
                    expression(paste(italic(' W'[peak]), sep = '')),
                    expression(paste(italic('    W'[eq]), sep = ''))), #col = 3,
         lty = 2,  mgp = c(3,1,0), cex.axis = 1.5)
    axis(1, at = c(-4,0,4), mgp = c(3,0.8,0), cex.axis = 1.6)
    axis(2, at = max(r0.seq), labels = expression(paste(italic('R'[0]), sep = '')), #col = 2, 
         lty=2, las = 2, mgp = c(3,0.75,0), cex.axis = 1.5)
    axis(2, at = c(-0.3,  -0.15, 0, 0.15,  0.3), cex.axis = 1.6, 
         labels = c('-0.3',  '-0.15', '0', '0.15',  '0.3'), mgp = c(3,0.5,0))
    
    mtext(side = 2, text = expression(italic(BBR)), line = 2.25)
    text(x=5.2, y=0.3, labels='a', pos=1, cex=3)
    
        
#plot dwdt trajectories ############
plot(log(W.seq), dwdt.seq, type = 'l', lwd = 2, ylim = c(-0.002, 0.01),
     xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
    
    mtext(side = 1, text = expression(paste('  log (', italic('W'), ')', sep = '')), line = 2.5)
    mtext(side = 2, text = expression(paste(italic(frac(dW, dt)))), line = 1.25)
    
    axis(1, at = c(-4,0,4), mgp = c(3,0.8,0), cex.axis = 1.6)
    axis(2, at = c(-0.002,0,0.002,0.004,0.006,0.008, 0.01),
         labels = c('-0.002','0','', '', '0.006','','0.01'), mgp = c(3,0.75,0), cex.axis = 1.6)
    
    lines(log(W.seq), dwdt.seq2, lwd = 2, lty = 2)
    text(x=5.2, y=0.01, labels='b', pos=1, cex=3)
    legend(x = 2, y = -0.001, legend = c(expression(italic(PDD-free)), expression(italic(PDD))),
           lwd = 2, lty = c(2,1), bty = 'n', cex = 1.2)
    
#zoom to y axis
plot(log(W.seq), dwdt.seq, type = 'l', lwd = 2, xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '',
     ylab = '', #xlab = expression(paste('log (', italic('W'), ')', sep = '')),
     ylim = c(-0.0002, 0.0002), xlim = c(-4,0))

    #mtext(side = 2, text = expression(paste(italic(frac(dW, dt)))), line = 1.5)
    axis(1, at = c(-4,-2,0), mgp = c(3,0.75,0), cex.axis = 1.6)
    axis(2, at = c(-0.0002, 0,0.0002),
         labels = c('-0.0002',  '0', '0.0002'),
         mgp = c(3,0.5,0), cex.axis = 1.6)
    
    
    segments(x0 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
             x1 = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]), #col = 3,
             y0 = 0, y1 = -1, lty = 2)
    
    axis(1, at = log(W.seq[which(round(Reff.seq, digits = 2) == 1.0)][1]),
         labels = expression(paste(italic('W'[bp]), sep = '')),
         lty = 2,  mgp = c(3,1,0), cex.axis = 1.5)
    
    abline(h = 0, lty=3)
    lines(log(W.seq), dwdt.seq2, lwd = 2, lty = 2)