#Toxicity to Biomphalaria snails from Ibrahim 1992 ###################
# reduction in snail fecundity: from table 1, relative change in number of juveniles produced over the 8 week experiment
snail.repro = data.frame(dose = c(0,125,250,500),
                         f_red = c(1, 3005/4225, 1749/4225, 1313/4225))


plot(snail.repro$dose, snail.repro$f_red, ylim = c(0,1), pch = 16, 
     ylab = 'reduction in snail recruitment rate', xlab = 'ChlorP ppb', main='Ibrahim92 snail recruitment')

mod1<-nls(f_red ~ exp(-b*(dose)), data=snail.repro, start = list(b=0.01))
  summary(mod1)

dose = c(0:500)

f_Nq_chlor_ibrahim92 = function(In){
  exp(-0.0028479*(In))
}

lines(dose, f_Nq_chlor_ibrahim92(dose), lty=2, col='red')
text(75,0.1, labels = expression(paste('f'[N],'(q) = ', 'e'^'-bq', '     b=0.00285', sep='')))

#Snail mortality -- relative change in survival over entire experiment period (12 weeks)
snail.mort = data.frame(dose = c(0,125,250,500),
                        mort = c(16.7, 23.3, 26.7, 100))

plot(snail.mort$dose, (snail.mort$mort - snail.mort$mort[1])/100, ylim = c(0,1), pch = 16, 
     ylab = 'increased snail mortality', xlab = 'ChlorP ppb', main = 'Ibrahim 92 snail mortality')

ibr_muNq<-drm((snail.mort$mort - snail.mort$mort[1])/100 ~ snail.mort$dose,
              data = snail.mort, type = 'binomial', fct = LL.2())

muNq_chlor_ibrahim92<-function(In){
  1/(1+exp(ibr_muNq$coefficients[1]*(log(In)-log(ibr_muNq$coefficients[2]))))
}  

lines(dose, muNq_chlor_ibrahim92(dose), lty=2, col='red')
text(150,0.5, labels = expression(paste(mu[N],'(q) = ', 'f(SF,LC'[50],',q)', '     SF=-3.87, LC'[50],'=358.5', sep='')), cex=0.7)
