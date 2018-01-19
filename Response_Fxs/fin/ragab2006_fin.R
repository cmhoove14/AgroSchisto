#Data extraction and model fitting to Ragab 2006 data
require(drc)

#Ammonium Nitrate Snail toxicity ##########
lc50.amm.rag = 470
slp.amm.rag = 10^1.30
se.lc50.amm.rag = mean(c(log10(507.6 / lc50.amm.rag), 
                         log10(lc50.amm.rag / 435.19))) / 1.96 
lc90.amm.rag = 640
  
fx.amm.rag = function(Fe, lc = lc50.amm.rag){
  feu = Fe/1000
  pnorm(slp.amm.rag * log10(feu / lc))
} 

plot(c(0,lc50.amm.rag, lc90.amm.rag)*1000, c(0,0.5, 0.9), pch = 16,
     xlim = c(0, (lc90.amm.rag+100)*1000), ylim = c(0,1),
     xlab = 'amm. phosphate (ppb)', ylab = 'mortality',
     main = expression(paste('Ragab06 Ammonium phosphate toxicity to ', 
                             italic('Bi. alexandrina'))))
  segments(x0 = 435.19*1000, x1 = 507.6*1000, y0 = 0.5, y1 = 0.5)

  lines(seq(0, (lc90.amm.rag+100)*1000, 750),
        sapply(seq(0, (lc90.amm.rag+100)*1000, 750), fx.amm.rag),
        lty = 2, col = 2)
  lines(seq(0, (lc90.amm.rag+100)*1000, 750),
        sapply(seq(0, (lc90.amm.rag+100)*1000, 750), fx.amm.rag, lc = 435.19),
        lty = 2, col = 2)
  lines(seq(0, (lc90.amm.rag+100)*1000, 750),
        sapply(seq(0, (lc90.amm.rag+100)*1000, 750), fx.amm.rag, lc = 507.60),
        lty = 2, col = 2)
  
rag06_mun_amm = function(Fe){
    feu = Fe/1000
    lc50 = 10^(rnorm(1, log10(lc50.amm.rag), se.lc50.amm.rag))
    mun = pnorm((slp.amm.rag) * log10(feu / lc50))
  
  return(mun)
}

points(seq(0, (lc90.amm.rag+100)*1000, 1000),
       sapply(seq(0, (lc90.amm.rag+100)*1000, 1000), rag06_mun_amm),
       pch = 5, col = 4, cex = 0.5)

#Pottasium Sulfate Snail toxicity ##########    
lc50.pot.rag = 1900
slp.pot.rag = 10^1.27
se.lc50.pot.rag = mean(c(log10(2280 / lc50.pot.rag), 
                         log10(lc50.pot.rag / 1583.3))) / 1.96 
lc90.pot.rag = 2600

fx.pot.rag = function(Fe, lc = lc50.pot.rag){
  feu = Fe/1000
  pnorm(slp.pot.rag * log10(feu / lc))
} 

plot(c(0, lc50.pot.rag, lc90.pot.rag)*1000, c(0, 0.5, 0.9), pch = 16,
     xlim = c(0, (lc90.pot.rag+100)*1000), ylim = c(0,1),
     xlab = 'pot. sulfate (ppb)', ylab = 'mortality',
     main = expression(paste('Ragab06 Potassium Sulfate toxicity to ', 
                             italic('Bi. alexandrina'))))
segments(x0 = 1583.3*1000, x1 = 2280*1000, y0 = 0.5, y1 = 0.5)

  lines(seq(0, (lc90.pot.rag+100)*1000, 750),
        sapply(seq(0, (lc90.pot.rag+100)*1000, 750), fx.pot.rag),
        lty = 2, col = 2)
  lines(seq(0, (lc90.pot.rag+100)*1000, 750),
        sapply(seq(0, (lc90.pot.rag+100)*1000, 750), fx.pot.rag, lc = 1583.3),
        lty = 2, col = 2)
  lines(seq(0, (lc90.pot.rag+100)*1000, 750),
        sapply(seq(0, (lc90.pot.rag+100)*1000, 750), fx.pot.rag, lc = 2280),
        lty = 2, col = 2)

rag06_mun_pot = function(Fe){
    feu = Fe/1000
    lc50 = 10^(rnorm(1, log10(lc50.pot.rag), se.lc50.pot.rag))
    mun = pnorm((slp.pot.rag) * log10(feu / lc50))
  
  return(mun)
}

points(seq(0, (lc90.pot.rag+100)*1000, 2000),
       sapply(seq(0, (lc90.pot.rag+100)*1000, 2000), rag06_mun_pot),
       pch = 5, col = 4, cex = 0.5)

#Urea Snail toxicity ##########    
lc50.urea.rag = 22000
slp.urea.rag = 10^1.28
se.lc50.urea.rag = mean(c(log10(24860 / lc50.urea.rag), 
                          log10(lc50.urea.rag / 19469))) / 1.96 
lc90.urea.rag = 31000

fx.urea.rag = function(Fe, lc = lc50.urea.rag){
  feu = Fe/1000
  pnorm(slp.urea.rag * log10(feu / lc))
} 

plot(c(0, lc50.urea.rag, lc90.urea.rag)*1000, c(0, 0.5, 0.9), pch = 16,
     xlim = c(0, (lc90.urea.rag+100)*1000), ylim = c(0,1),
     xlab = 'urea (ppb)', ylab = 'mortality',
     main = expression(paste('Ragab06 urea toxicity to ', 
                             italic('Bi. alexandrina'))))
  segments(x0 = 24860*1000, x1 = 19469*1000, y0 = 0.5, y1 = 0.5)

  lines(seq(0, (lc90.urea.rag+100)*1000, 750),
        sapply(seq(0, (lc90.urea.rag+100)*1000, 750), fx.urea.rag),
        lty = 2, col = 2)
  lines(seq(0, (lc90.urea.rag+100)*1000, 750),
        sapply(seq(0, (lc90.urea.rag+100)*1000, 750), fx.urea.rag, lc = 19469),
        lty = 2, col = 2)
  lines(seq(0, (lc90.urea.rag+100)*1000, 750),
        sapply(seq(0, (lc90.urea.rag+100)*1000, 750), fx.urea.rag, lc = 24860),
        lty = 2, col = 2)

rag06_mun_urea = function(Fe){
    feu = Fe/1000
    lc50 = 10^(rnorm(1, log10(lc50.urea.rag), se.lc50.urea.rag))
    mun = pnorm((slp.urea.rag) * log10(feu / lc50))
  
  return(mun)
}

points(seq(0, (lc90.urea.rag+100)*1000, 20000),
       sapply(seq(0, (lc90.urea.rag+100)*1000, 20000), rag06_mun_urea),
       pch = 5, col = 4, cex = 0.5)

#Keep vector #################
keep.ragab.mun = c('rag06_mun_urea', 'lc50.urea.rag', 'se.lc50.urea.rag', 'slp.urea.rag',
                   'rag06_mun_pot', 'lc50.pot.rag', 'se.lc50.pot.rag', 'slp.pot.rag',
                   'rag06_mun_amm', 'lc50.amm.rag', 'se.lc50.amm.rag', 'slp.amm.rag')
