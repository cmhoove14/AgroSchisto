#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#load workspace
load('Review_models/r0_of_malatr_sim.RData') #conc range = 0 - 2000 by 2s, 2000 - 2e5 by 200s

#pre-process of the dataframe
  df.malatra$atr.norm = df.malatra$atr / max(df.malatra$atr)
  df.malatra$mal.norm = df.malatra$mal / max(df.malatra$mal)

#3d plot of mean change in R0 across all sims
fin.malatr = matrix(df.malatra[,(nsims+3)], ncol = nconc)

rg = colorRampPalette(c('green', 'red'))
  colmalatr = rg(100)
  zfacet = (fin.malatr[-1, -1] + 
              fin.malatr[-1, -ncol(fin.malatr)] + 
              fin.malatr[-nrow(fin.malatr), -1] + 
              fin.malatr[-nrow(fin.malatr), -ncol(fin.malatr)]) / 4
  facetcol = cut(zfacet, 100)

persp(y = mal.range, ylim = range(mal.range), x = atr.range, xlim = range(atr.range),
      z = fin.malatr - r0.Ag(0)[3], ticktype = 'simple', nticks = 4, zlim = c(-4.1,4.1),
      ylab = 'Malathion', xlab = 'Atrazine', zlab = 'Change in R0',
      phi = 20, theta = 45, col = colmalatr[facetcol])

#Plot smoothed version  
mod3d = gam(mean ~ te(atr, mal, k = 20), data = df.malatra)  

  smth.malatr = matrix(fitted(mod3d), ncol = nconc)

  zfacet2 = (smth.malatr[-1, -1] + 
              smth.malatr[-1, -ncol(smth.malatr)] + 
              smth.malatr[-nrow(smth.malatr), -1] + 
              smth.malatr[-nrow(smth.malatr), -ncol(smth.malatr)]) / 4
  facetcol2 = cut(zfacet2, 100)

persp(y = mal.range, ylim = range(mal.range), x = atr.range, xlim = range(atr.range),
      z = smth.malatr - r0.Ag(0)[3], ticktype = 'detailed', nticks = 4, zlim = c(-4.1,4.1),
      ylab = 'Malathion', xlab = 'Atrazine', zlab = 'Change in R0',
      phi = 20, theta = 60, col = colmalatr[facetcol2])

#Find volume of concentration space resulting in increased R0
  df.malatra.d = df.malatra
    df.malatra.d[,c(3:(nsims+2))] = df.malatra.d[,c(3:(nsims+2))] - r0.Ag(0)[3]
    
  df.malatra.d0 = df.malatra.d
    df.malatra.d0[df.malatra.d0 < 0] = 0
    
  da = 1/(nconc-1)
    a = unique(df.malatra.d$atr.norm)
  dm = 1/(nconc-1)
    m = unique(df.malatra.d$mal.norm)
    
  vr0 = array(data = NA, dim = c((nconc - 1), (nconc - 1), nsims, 2))

  for(h in 3:(nsims+2)){
   for(i in 1:(nconc - 1)){
    for(j in 1:(nconc - 1)){
      print(c(i,j,h))
      
      h1 = mean(df.malatra.d[df.malatra.d$atr.norm == a[i] & df.malatra.d$mal.norm == m[j],h],
                df.malatra.d[df.malatra.d$atr.norm == a[i] & df.malatra.d$mal.norm == m[j+1],h],
                df.malatra.d[df.malatra.d$atr.norm == a[i+1] & df.malatra.d$mal.norm == m[j],h],
                df.malatra.d[df.malatra.d$atr.norm == a[i+1] & df.malatra.d$mal.norm == m[j+1],h])
      
      h2 = mean(df.malatra.d0[df.malatra.d0$atr.norm == a[i] & df.malatra.d0$mal.norm == m[j],h],
                df.malatra.d0[df.malatra.d0$atr.norm == a[i] & df.malatra.d0$mal.norm == m[j+1],h],
                df.malatra.d0[df.malatra.d0$atr.norm == a[i+1] & df.malatra.d0$mal.norm == m[j],h],
                df.malatra.d0[df.malatra.d0$atr.norm == a[i+1] & df.malatra.d0$mal.norm == m[j+1],h])
      
      fill1 = da*dm*h1
      fill2 = da*dm*h2
      vr0[i,j,h-2,1] = fill1
      vr0[i,j,h-2,2] = fill2
    }
  }
}  

#Volume of entire surface   
  vs = as.numeric()
  for(i in 1:nsims){
    vs[i] = sum(vr0[, , i, 1])
  }

#Volume of increased r0 area    
  vs0 = as.numeric()
  for(i in 1:nsims){
    vs0[i] = sum(vr0[, , i, 2])
  }
 
#Smoothed area of r0>0 ############
  smth.malatrg0 = smth.malatr - r0.Ag(0)[3] 
  smth.malatrg0[smth.malatrg0 < 0] = 0
  
  zfacet3 = (smth.malatrg0[-1, -1] + 
               smth.malatrg0[-1, -ncol(smth.malatrg0)] + 
               smth.malatrg0[-nrow(smth.malatrg0), -1] + 
               smth.malatrg0[-nrow(smth.malatrg0), -ncol(smth.malatrg0)]) / 4
  facetcol3 = cut(zfacet3, 100)
  
  persp(y = mal.range, ylim = range(mal.range), x = atr.range, xlim = range(atr.range),
        z = smth.malatrg0, ticktype = 'detailed', nticks = 4, zlim = c(0,2),
        ylab = 'Malathion', xlab = 'Atrazine', zlab = 'Change in R0',
        phi = 30, theta = 60, col = colmalatr[facetcol3])

