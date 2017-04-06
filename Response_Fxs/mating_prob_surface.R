k.vec = seq(1e-10, 2, length.out = 100)
w.vec = seq(1e-12,80, length.out = 100)

fx<-function(x, w, k){
  (1-cos(x)) / ((1+((w/(w+k))*cos(x)))^(1+k))
}

get.phi = function(w,k){
  alph = w / (w+k)
  p1 = ((1-alph)^(1+k)) / (2*pi)
  p2 = integrate(fx, 0, (2*pi), w, k, stop.on.error = F)$value
  phi = 1 - p1 * p2
  
  return(c(p1, p2, phi))
}

mat.init = array(data = NA, dim = c(100,100,3))

for(i in 1:100){
  mat.init[ , i, 1] = sapply(k.vec, get.phi, w = w.vec[i], simplify = TRUE)[1,]
  mat.init[ , i, 2] = sapply(k.vec, get.phi, w = w.vec[i], simplify = TRUE)[2,]
  mat.init[ , i, 3] = sapply(k.vec, get.phi, w = w.vec[i], simplify = TRUE)[3,]
  
  print(i)
}

mat.p1 = mat.init[ , , 1]
mat.p2 = mat.init[ -1, -1, 2]
mat.phi = mat.init[ , , 3]

windows(record = T)
persp(y = w.vec, ylim = range(w.vec), x = k.vec, xlim = range(k.vec),
      z = mat.phi, zlim = c(0,1.1), ticktype = 'detailed', nticks = 4, 
      ylab = 'mean worm burden',
      xlab = 'clumping parameter',
      zlab = 'mating probability',
      phi = 20, theta = 45, shade = 0.4, col = 'lightblue')

persp(y = w.vec, ylim = range(w.vec), x = k.vec, xlim = range(k.vec),
      z = mat.p1, ticktype = 'detailed', nticks = 4, 
      ylab = 'mean worm burden',
      xlab = 'clumping parameter',
      zlab = 'p1',
      phi = 20, theta = 45, shade = 0.4, col = 'lightblue')

persp(y = w.vec[-1], ylim = range(w.vec[-1]), x = k.vec[-1], xlim = range(k.vec[-1]),
      z = mat.p2, ticktype = 'detailed', nticks = 4, 
      ylab = 'mean worm burden',
      xlab = 'clumping parameter',
      zlab = 'p2',
      phi = 20, theta = 45, shade = 0.4, col = 'lightblue')

#Fit simplified function to mating function for use in stochastic modeling code ########
phi.parms = matrix(ncol = 3, nrow = length(k.vec))
  phi.parms[,1] = k.vec
  
plot(w.vec, mat.phi[1,], type = 'l', lwd = 2)
  for(i in 2:length(k.vec)){
    lines(w.vec, mat.phi[i,], col = i)
  }

for(i in 1:length(k.vec)){
  phi.parms[i,c(2,3)] = summary(nls(mat.phi[i,] ~ a-a/(w.vec+1)^(b), start = list(a=1,b = 1),
                                    upper = list(a=1,b = Inf), 
                                    lower = list(a=0,b = 0), 
                                    algorithm = 'port'))$coefficients[1:2]
}

phi.simp = function(a,b,w){
  a-a/(w+1)^(b)
}

plot(w.vec, sapply(w.vec, phi.simp, a = phi.parms[1,2], b = phi.parms[1,3], simplify = T),
     type = 'l', lwd = 2)
for(i in 2:length(k.vec)){
  lines(w.vec, sapply(w.vec, phi.simp, a = phi.parms[i,2], b = phi.parms[i,3], simplify = T), col = i)
}

colnames(phi.parms) = c('k', 'a', 'b')
write.csv(phi.parms, 'Elimination_Feasibility/phi_parameters_simplified.csv', row.names = F)