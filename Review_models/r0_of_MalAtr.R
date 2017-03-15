#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Load response and R0 functions #######
source('Response_Fxs/bakry2012.R')
source('Response_Fxs/piC_atr_meta_beq.R')
source('Response_Fxs/koprivnikar06_piC_beq.R')
source('Response_Fxs/griggs08_piC_beq.R')
source('Response_Fxs/Baxter_Rohr_Atrazine2011.R')
source('Response_Fxs/rohr08_piC_atrOnly.R')
source('Response_Fxs/tchounwou92_piC2_beq.R')
source('Response_Fxs/tchounwou91_piM_beq.R')
source('Response_Fxs/Halstead_Insecticides2015.R')
source('Response_Fxs/bakry2011.R')
source('Response_Fxs/tchounwou91_fN-muN.R')

source('Review_models/r0_of_q.R')

library(parallel)
library(reshape2)
library(mgcv)

keep.fin.malatr = c(keep.bak12, keep.grg08, keep.kop06.beq, keep.meta.piC, keep.baxrohr, keep.atr.rohr08,
                 'r0.He', 'r0.fix', 'r0.In', 'r0.Ag', 'keep.fin.malatr', 'parameters', 'nil0', 'nil1', 
                 keep.tch92.beq, keep.tch91.beq, keep.tch91.snail, keep.hal15.muP[c(1,7,13)], 
                 keep.bak11.N[c(1,2,3,6,7)])

rm(list = setdiff(ls(), keep.fin.malatr))
dev.off()

no.cores = detectCores() - 1

#3d plot of atrazine/malathion mixture on R0 #################
nsims = 30  
nconc = 100
rec.mal = 100000
  mal.range = seq(0, rec.mal, length.out = nconc)
  mal.vec = rep(mal.range, each = nconc)
rec.atr = 300
  atr.range = seq(0, rec.atr, length.out = nconc)
  atr.vec = rep(atr.range, times = nconc)

m.malatr = matrix(ncol = nsims+2, nrow = nconc^2) 
  m.malatr[,1] = atr.vec
  m.malatr[,2] = mal.vec

clusmalatr = makeCluster(no.cores)
  clusterExport(clusmalatr, c(keep.fin.malatr, 'uniroot.all', 'rdrm', 'LL.2'))  
  
  for(i in 1:nsims){
    set.seed(i)
    m.malatr[,(i+2)] = clusterMap(clusmalatr, fun = r0.Ag, He = atr.vec, In = mal.vec, SIMPLIFY = TRUE,
                                  MoreArgs = c(f.in.f_Nq = fNq_mal_tch91_uncertainty, 
                                               f.in.mu_Pq = muPq_mal_Halstead_uncertainty,
                                               f.he.phi_Nq = phi_Nq_atr_baxrohr.no30, 
                                               f.in.mu_Nq = muNq_mal_tch91_uncertainty, 
                                               f.he.mu_Nq = muNq_atr_bak12, 
                                               f.in.pi_Mq = piM.tch91_mal_unc, 
                                               f.in.pi_Cq = piC.tch92_mal_unc, 
                                               f.he.pi_Cq = piC.atr.rohr08.lin))[3,]
  }

stopCluster(clusmalatr)

#Bit of post-processing ##########  
  df.malatra = as.data.frame(m.malatr)

colnames(df.malatra)[1:2] = c('atr', 'mal')
  
  df.malatra$mean = rowMeans(df.malatra[,3:nsims])
  
#Find ~volume of concentration space resulting in increased R0
  df.malatra.d = df.malatra
  df.malatra.d[,c(3:(nsims+2))] = df.malatra.d[,c(3:(nsims+2))] - r0.Ag(0)[3]
  
  df.malatra.d0 = df.malatra.d
  df.malatra.d0[df.malatra.d0 < 0] = 0
  
  da = 1/(nconc-1) #change in atrazine concentration = length
    a = unique(df.malatra.d$atr.norm)
  dm = 1/(nconc-1) #change in malathion concentration = width
    m = unique(df.malatra.d$mal.norm)
  
  vr0 = array(data = NA, dim = c((nconc - 1), (nconc - 1), nsims, 2))
  
  for(h in 3:(nsims+2)){
    for(i in 1:(nconc - 1)){
      for(j in 1:(nconc - 1)){
      #mean of change in r0 at four concentration cooredinates approximates height
        h1 = mean(df.malatra.d[df.malatra.d$atr.norm == a[i] & df.malatra.d$mal.norm == m[j],h],
                  df.malatra.d[df.malatra.d$atr.norm == a[i] & df.malatra.d$mal.norm == m[j+1],h],
                  df.malatra.d[df.malatra.d$atr.norm == a[i+1] & df.malatra.d$mal.norm == m[j],h],
                  df.malatra.d[df.malatra.d$atr.norm == a[i+1] & df.malatra.d$mal.norm == m[j+1],h])
        
        h2 = mean(df.malatra.d0[df.malatra.d0$atr.norm == a[i] & df.malatra.d0$mal.norm == m[j],h],
                  df.malatra.d0[df.malatra.d0$atr.norm == a[i] & df.malatra.d0$mal.norm == m[j+1],h],
                  df.malatra.d0[df.malatra.d0$atr.norm == a[i+1] & df.malatra.d0$mal.norm == m[j],h],
                  df.malatra.d0[df.malatra.d0$atr.norm == a[i+1] & df.malatra.d0$mal.norm == m[j+1],h])
        
        fill1 = da*dm*h1 #V = length*width*height
        fill2 = da*dm*h2 #V = length*width*height
        vr0[i,j,h-2,1] = fill1
        vr0[i,j,h-2,2] = fill2
      }
    }
  }  
  