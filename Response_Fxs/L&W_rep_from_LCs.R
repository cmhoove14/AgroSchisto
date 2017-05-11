require(LW1949)
require(drc)

conc <- c(0.0625, 0.125, 0.25, 0.5, 1, 2, 3)
numtested <- rep(8, 7)
numaffected <- c(1, 4, 4, 7, 8, 8, 8)
mydat <- dataprep(dose=conc, ntot=numtested, nfx=numaffected)

  plot(mydat$dose, mydat$pfx, pch = 16, xlim = c(0,3), ylim = c(0,1),
       xlab = 'Dose', ylab = 'mortality')
  
intslope <- fitLWauto(mydat)

plot(mydat$log10dose, mydat$bitpfx, pch = 16, xlab = 'log10 dose', ylab = 'prob(mortality)')
  abline(intslope, lty = 2, col = 2)
  
fLW <- LWestimate(intslope, mydat)

lcs = as.data.frame(predlinear(c(10,25,50,90), fLW))
  lcs$probit = qnorm(lcs$pct/100)
  lcs$log10 = log10(lcs$ED)
plot(lcs$log10, lcs$probit, pch = 16, ylim = c(-1.5,1.5),
     xlim = c(-1.75,0.5), xlab = 'log10(dose)', ylab = 'probit(mortality)')
  lm.lcs = lm(probit ~ log10, data = lcs)
  abline(coef(lm.lcs), lty = 2, col = 2)
#able to recover slope and intercept from just two points?  
    lcs2 = subset(lcs, probit >= 0)
    lm.lcs2 = lm(probit ~ log10, data = lcs2)
  abline(coef(lm.lcs2), lty = 2, col = 3, lwd = 2)  
#yes
    predLinesLP(fLW)
  segments(x0 = log10(lcs[,3]), x1 = log10(lcs[,4]), 
           y0 = qnorm(lcs[,1]/100), y1 = qnorm(lcs[,1]/100))  
  
plot(mydat$dose, mydat$pfx, pch = 16, xlim = c(0,3), ylim = c(0,1),
     xlab = 'Dose', ylab = 'mortality')
  lines(predlinear(c(1:99), fLW)[,2], predlinear(c(1:99), fLW)[,1]/100, lty = 2, col = 2)
  lines(predlinear(c(1:99), fLW)[,3], predlinear(c(1:99), fLW)[,1]/100, lty = 3, col = 2)
  lines(predlinear(c(1:99), fLW)[,4], predlinear(c(1:99), fLW)[,1]/100, lty = 3, col = 2)
  
fx = function(d, lc){
  pnorm(coef(lm.lcs)[2] * log10(d / lc))
}  

  lines(seq(0,3,0.01), sapply(seq(0,3,0.01), fx, lc = lcs$ED[3]), lty = 2)
  lines(seq(0,3,0.01), sapply(seq(0,3,0.01), fx, lc = lcs$lower[3]), lty = 2)
  lines(seq(0,3,0.01), sapply(seq(0,3,0.01), fx, lc = lcs$upper[3]), lty = 2)
#doesn't exactly get uncertainty right, but it's fairly close
  
#slope function by hand
  slp = predlinear(c(16,50,84), fLW)