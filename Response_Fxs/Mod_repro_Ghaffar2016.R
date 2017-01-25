require(LW1949)
require(drc)

lc50.n.but = 5.560
slp.n.but = 1.093

lcs = c(10, 25, 50, 90)/100
  probs = probit(lcs/100)
lc.but.vals = c(2417, 3906, 5560, 8703)
  concs = log10(lc.but.vals/1000)
    
plot(concs, probs, pch = 16, xlim = c(-1,1), ylim = c(-6,0))
  lm(probs ~ concs)
dose <- lc.but.vals/1000
ntested <- rep(10, 4)
nalive <- c(10, 25, 50, 90)/10
mydat <- dataprep(dose=dose, ntot=ntested, nfx=nalive)
fLW <- LWestimate(fitLWauto(mydat), mydat)
predlinear(c(10, 25, 50, 90), fLW)[1:4,2]

plot(10^concs, lcs*100, pch = 17, ylim = c(0,100), xlim = c(0,10))

predLines(fLW)

points(predlinear(c(10, 25, 50, 90), fLW)[1:4,2], predlinear(c(10, 25, 50, 90), fLW)[1:4,1],  pch = 16, col = 'red')

predlinear(c(10, 25, 50, 90), c(-3.160386, 4.469945), simple = TRUE)

lines(predlinear(c(1:99), c(-3.160386, 4.469945), simple = TRUE), c(1:99),lty=2, col=2)