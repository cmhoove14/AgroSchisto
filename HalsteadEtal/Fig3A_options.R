#Step 2.1.1 Get R0 estimates across agrochemical treatment groups for low transmission#####################
#village R0        
R0_vil.low = get_Ro_beta_lamda(beta = pt2$beta, lamda = pt2$lamda, muPq = p.dead)[3]
R0_up.low = get_Ro_beta_lamda(beta = pt2$beta, lamda = pt2$lamda, muPq = p.dead)[3]
R0_lo.low = get_Ro_beta_lamda(beta = pt2$beta, lamda = pt2$lamda, muPq = p.dead)[3]

#village R0 with prawns
R0_vil.p.low<-get_Ro_beta_lamda(beta=pt2$beta, lamda = pt2$lamda)[3]
R0_up.p.low<-get_Ro_beta_lamda(beta=pt2$beta, lamda = pt2$lamda)[3]
R0_lo.p.low<-get_Ro_beta_lamda(beta=pt2$beta, lamda = pt2$lamda)[3]


#atrazine only R0
R0_atra0.low = get_Ro_beta_lamda(phi_Nq = phi_Nqs[3], 
                                 #f_Nq = f_Nqs[3], 
                                 #f_Nq = fNqs.ch[3],
                                 beta = pt2$beta, 
                                 lamda = pt2$lamda)[3]
R0_atra_up.low = get_Ro_beta_lamda(phi_Nq = phi_Nqs.hi[3], 
                                   #f_Nq = f_Nqs[3], 
                                   #f_Nq = fNqs.ch[3],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3]
R0_atra_lo.low = get_Ro_beta_lamda(phi_Nq = phi_Nqs.lo[3], 
                                   #f_Nq = f_Nqs[3], 
                                   #f_Nq = fNqs.ch[3],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3]

#chlorpyrifos only R0
R0_chlor0.low = get_Ro_beta_lamda(muPq = muPq, 
                                  beta = pt2$beta, 
                                  lamda = pt2$lamda,
                                  phi_Nq = phi_Nqs[1])[3]
R0_chlor_up.low = get_Ro_beta_lamda(muPq = muPq, 
                                    beta = pt2$beta, 
                                    lamda = pt2$lamda,
                                    phi_Nq = phi_Nqs.hi[1])[3]
R0_chlor_lo.low = get_Ro_beta_lamda(muPq = muPq, 
                                    beta = pt2$beta, 
                                    lamda = pt2$lamda,
                                    phi_Nq = phi_Nqs.lo[1])[3]

#fertilizer only R0
R0_fert0.low = get_Ro_beta_lamda(phi_Nq = phi_Nqs[2], 
                                 #f_Nq = f_Nqs[2], 
                                 #f_Nq = fNqs.ch[2],
                                 beta = pt2$beta, 
                                 lamda = pt2$lamda)[3]
R0_fert_up.low = get_Ro_beta_lamda(phi_Nq = phi_Nqs.hi[2], 
                                   #f_Nq = f_Nqs[2], 
                                   #f_Nq = fNqs.ch[2],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3]
R0_fert_lo.low = get_Ro_beta_lamda(phi_Nq = phi_Nqs.lo[2], 
                                   #f_Nq = f_Nqs[2],
                                   #f_Nq = fNqs.ch[2],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3]  

#atrazine+chlorpyrifos R0
R0_atch0.low = get_Ro_beta_lamda(muPq = muPq, 
                                 phi_Nq = phi_Nqs[3], 
                                 #f_Nq = f_Nqs[3], 
                                 #f_Nq = fNqs.ch[3],
                                 beta = pt2$beta, 
                                 lamda = pt2$lamda)[3]
R0_atch_up.low = get_Ro_beta_lamda(muPq = muPq, 
                                   phi_Nq = phi_Nqs.hi[3], 
                                   #f_Nq = f_Nqs[3], 
                                   #f_Nq = fNqs.ch[3],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3]
R0_atch_lo.low = get_Ro_beta_lamda(muPq = muPq, 
                                   phi_Nq = phi_Nqs.lo[3], 
                                   #f_Nq = f_Nqs[3], 
                                   #f_Nq = fNqs.ch[3],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3] 

#atrazine+fertilizer R0
R0_atfe0.low = get_Ro_beta_lamda(phi_Nq = phi_Nqs[4], 
                                 #f_Nq = f_Nqs[4], 
                                 #f_Nq = fNqs.ch[4],
                                 beta = pt2$beta, 
                                 lamda = pt2$lamda)[3]
R0_atfe_up.low = get_Ro_beta_lamda(phi_Nq = phi_Nqs.hi[4], 
                                   #f_Nq = f_Nqs[4], 
                                   #f_Nq = fNqs.ch[4],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3]
R0_atfe_lo.low = get_Ro_beta_lamda(phi_Nq = phi_Nqs.lo[4], 
                                   #f_Nq = f_Nqs[4], 
                                   #f_Nq = fNqs.ch[4],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3] 

#chlorpyrifos+fertilizer R0
R0_chfe0.low = get_Ro_beta_lamda(muPq = muPq, 
                                 phi_Nq = phi_Nqs[2], 
                                 #f_Nq = f_Nqs[2],
                                 #f_Nq = fNqs.ch[2],
                                 beta = pt2$beta, 
                                 lamda = pt2$lamda)[3]
R0_chfe_up.low = get_Ro_beta_lamda(muPq = muPq, 
                                   phi_Nq = phi_Nqs.hi[2], 
                                   #f_Nq = f_Nqs[2], 
                                   #f_Nq = fNqs.ch[2],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3]
R0_chfe_lo.low = get_Ro_beta_lamda(muPq = muPq, 
                                   phi_Nq = phi_Nqs.lo[2], 
                                   #f_Nq = f_Nqs[2], 
                                   #f_Nq = fNqs.ch[2],
                                   beta = pt2$beta, 
                                   lamda = pt2$lamda)[3]  

#atrazine+chlorpyrifos+fertilizer R0
R0_tre0.low = get_Ro_beta_lamda(muPq = muPq, 
                                phi_Nq = phi_Nqs[4], 
                                #f_Nq = f_Nqs[4], 
                                #f_Nq = fNqs.ch[4],
                                beta = pt2$beta, 
                                lamda = pt2$lamda)[3]
R0_tre_up.low = get_Ro_beta_lamda(muPq = muPq, 
                                  phi_Nq = phi_Nqs.hi[4], 
                                  #f_Nq = f_Nqs[4], 
                                  #f_Nq = fNqs.ch[4],
                                  beta = pt2$beta, 
                                  lamda = pt2$lamda)[3]
R0_tre_lo.low = get_Ro_beta_lamda(muPq = muPq, 
                                  phi_Nq = phi_Nqs.lo[4], 
                                  #f_Nq = f_Nqs[4], 
                                  #f_Nq = fNqs.ch[4],
                                  beta = pt2$beta, 
                                  lamda = pt2$lamda)[3] 

#Step 2.1.2 Get R0 estimates across agrochemical treatment groups for high transmission#####################
#village R0        
R0_vil.mid = get_Ro_beta_lamda(beta = pt3$beta, lamda = pt3$lamda, muPq = p.dead)[3]
R0_up.mid = get_Ro_beta_lamda(beta = pt3$beta, lamda = pt3$lamda, muPq = p.dead)[3]
R0_lo.mid = get_Ro_beta_lamda(beta = pt3$beta, lamda = pt3$lamda, muPq = p.dead)[3]

#village R0 with prawns
R0_vil.p.mid<-get_Ro_beta_lamda(beta=pt3$beta, lamda = pt3$lamda)[3]
R0_up.p.mid<-get_Ro_beta_lamda(beta=pt3$beta, lamda = pt3$lamda)[3]
R0_lo.p.mid<-get_Ro_beta_lamda(beta=pt3$beta, lamda = pt3$lamda)[3]


#atrazine only R0
R0_atra0.mid = get_Ro_beta_lamda(phi_Nq = phi_Nqs[3], 
                                 #f_Nq = f_Nqs[3], 
                                 #f_Nq = fNqs.ch[3],
                                 beta = pt3$beta, 
                                 lamda = pt3$lamda)[3]
R0_atra_up.mid = get_Ro_beta_lamda(phi_Nq = phi_Nqs.hi[3], 
                                   #f_Nq = f_Nqs[3], 
                                   #f_Nq = fNqs.ch[3],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3]
R0_atra_lo.mid = get_Ro_beta_lamda(phi_Nq = phi_Nqs.lo[3], 
                                   #f_Nq = f_Nqs[3], 
                                   #f_Nq = fNqs.ch[3],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3]

#chlorpyrifos only R0
R0_chlor0.mid = get_Ro_beta_lamda(muPq = muPq, 
                                  beta = pt3$beta, 
                                  lamda = pt3$lamda,
                                  phi_Nq = phi_Nqs[1])[3]
R0_chlor_up.mid = get_Ro_beta_lamda(muPq = muPq, 
                                    beta = pt3$beta, 
                                    lamda = pt3$lamda,
                                    phi_Nq = phi_Nqs.hi[1])[3]
R0_chlor_lo.mid = get_Ro_beta_lamda(muPq = muPq, 
                                    beta = pt3$beta, 
                                    lamda = pt3$lamda,
                                    phi_Nq = phi_Nqs.lo[1])[3]

#fertilizer only R0
R0_fert0.mid = get_Ro_beta_lamda(phi_Nq = phi_Nqs[2], 
                                 #f_Nq = f_Nqs[2], 
                                 #f_Nq = fNqs.ch[2],
                                 beta = pt3$beta, 
                                 lamda = pt3$lamda)[3]
R0_fert_up.mid = get_Ro_beta_lamda(phi_Nq = phi_Nqs.hi[2], 
                                   #f_Nq = f_Nqs[2], 
                                   #f_Nq = fNqs.ch[2],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3]
R0_fert_lo.mid = get_Ro_beta_lamda(phi_Nq = phi_Nqs.lo[2], 
                                   #f_Nq = f_Nqs[2],
                                   #f_Nq = fNqs.ch[2],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3]  

#atrazine+chlorpyrifos R0
R0_atch0.mid = get_Ro_beta_lamda(muPq = muPq, 
                                 phi_Nq = phi_Nqs[3], 
                                 #f_Nq = f_Nqs[3], 
                                 #f_Nq = fNqs.ch[3],
                                 beta = pt3$beta, 
                                 lamda = pt3$lamda)[3]
R0_atch_up.mid = get_Ro_beta_lamda(muPq = muPq, 
                                   phi_Nq = phi_Nqs.hi[3], 
                                   #f_Nq = f_Nqs[3], 
                                   #f_Nq = fNqs.ch[3],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3]
R0_atch_lo.mid = get_Ro_beta_lamda(muPq = muPq, 
                                   phi_Nq = phi_Nqs.lo[3], 
                                   #f_Nq = f_Nqs[3], 
                                   #f_Nq = fNqs.ch[3],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3] 

#atrazine+fertilizer R0
R0_atfe0.mid = get_Ro_beta_lamda(phi_Nq = phi_Nqs[4], 
                                 #f_Nq = f_Nqs[4], 
                                 #f_Nq = fNqs.ch[4],
                                 beta = pt3$beta, 
                                 lamda = pt3$lamda)[3]
R0_atfe_up.mid = get_Ro_beta_lamda(phi_Nq = phi_Nqs.hi[4], 
                                   #f_Nq = f_Nqs[4], 
                                   #f_Nq = fNqs.ch[4],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3]
R0_atfe_lo.mid = get_Ro_beta_lamda(phi_Nq = phi_Nqs.lo[4], 
                                   #f_Nq = f_Nqs[4], 
                                   #f_Nq = fNqs.ch[4],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3] 

#chlorpyrifos+fertilizer R0
R0_chfe0.mid = get_Ro_beta_lamda(muPq = muPq, 
                                 phi_Nq = phi_Nqs[2], 
                                 #f_Nq = f_Nqs[2],
                                 #f_Nq = fNqs.ch[2],
                                 beta = pt3$beta, 
                                 lamda = pt3$lamda)[3]
R0_chfe_up.mid = get_Ro_beta_lamda(muPq = muPq, 
                                   phi_Nq = phi_Nqs.hi[2], 
                                   #f_Nq = f_Nqs[2], 
                                   #f_Nq = fNqs.ch[2],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3]
R0_chfe_lo.mid = get_Ro_beta_lamda(muPq = muPq, 
                                   phi_Nq = phi_Nqs.lo[2], 
                                   #f_Nq = f_Nqs[2], 
                                   #f_Nq = fNqs.ch[2],
                                   beta = pt3$beta, 
                                   lamda = pt3$lamda)[3]  

#atrazine+chlorpyrifos+fertilizer R0
R0_tre0.mid = get_Ro_beta_lamda(muPq = muPq, 
                                phi_Nq = phi_Nqs[4], 
                                #f_Nq = f_Nqs[4], 
                                #f_Nq = fNqs.ch[4],
                                beta = pt3$beta, 
                                lamda = pt3$lamda)[3]
R0_tre_up.mid = get_Ro_beta_lamda(muPq = muPq, 
                                  phi_Nq = phi_Nqs.hi[4], 
                                  #f_Nq = f_Nqs[4], 
                                  #f_Nq = fNqs.ch[4],
                                  beta = pt3$beta, 
                                  lamda = pt3$lamda)[3]
R0_tre_lo.mid = get_Ro_beta_lamda(muPq = muPq, 
                                  phi_Nq = phi_Nqs.lo[4], 
                                  #f_Nq = f_Nqs[4], 
                                  #f_Nq = fNqs.ch[4],
                                  beta = pt3$beta, 
                                  lamda = pt3$lamda)[3] 

#Step 2.1.3 Get R0 estimates across agrochemical treatment groups for intermediate transmission#####################
#village R0        
R0_vil.hi = get_Ro_beta_lamda(beta = pt4$beta, lamda = pt4$lamda, muPq = p.dead)[3]
R0_up.hi = get_Ro_beta_lamda(beta = pt4$beta, lamda = pt4$lamda, muPq = p.dead)[3]
R0_lo.hi = get_Ro_beta_lamda(beta = pt4$beta, lamda = pt4$lamda, muPq = p.dead)[3]

#village R0 with prawns
R0_vil.p.hi<-get_Ro_beta_lamda(beta=pt4$beta, lamda = pt4$lamda)[3]
R0_up.p.hi<-get_Ro_beta_lamda(beta=pt4$beta, lamda = pt4$lamda)[3]
R0_lo.p.hi<-get_Ro_beta_lamda(beta=pt4$beta, lamda = pt4$lamda)[3]


#atrazine only R0
R0_atra0.hi = get_Ro_beta_lamda(phi_Nq = phi_Nqs[3], 
                                #f_Nq = f_Nqs[3], 
                                #f_Nq = fNqs.ch[3],
                                beta = pt4$beta, 
                                lamda = pt4$lamda)[3]
R0_atra_up.hi = get_Ro_beta_lamda(phi_Nq = phi_Nqs.hi[3], 
                                  #f_Nq = f_Nqs[3], 
                                  #f_Nq = fNqs.ch[3],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3]
R0_atra_lo.hi = get_Ro_beta_lamda(phi_Nq = phi_Nqs.lo[3], 
                                  #f_Nq = f_Nqs[3], 
                                  #f_Nq = fNqs.ch[3],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3]

#chlorpyrifos only R0
R0_chlor0.hi = get_Ro_beta_lamda(muPq = muPq, 
                                 beta = pt4$beta, 
                                 lamda = pt4$lamda,
                                 phi_Nq = phi_Nqs[1])[3]
R0_chlor_up.hi = get_Ro_beta_lamda(muPq = muPq, 
                                   beta = pt4$beta, 
                                   lamda = pt4$lamda,
                                   phi_Nq = phi_Nqs.hi[1])[3]
R0_chlor_lo.hi = get_Ro_beta_lamda(muPq = muPq, 
                                   beta = pt4$beta, 
                                   lamda = pt4$lamda,
                                   phi_Nq = phi_Nqs.lo[1])[3]

#fertilizer only R0
R0_fert0.hi = get_Ro_beta_lamda(phi_Nq = phi_Nqs[2], 
                                #f_Nq = f_Nqs[2], 
                                #f_Nq = fNqs.ch[2],
                                beta = pt4$beta, 
                                lamda = pt4$lamda)[3]
R0_fert_up.hi = get_Ro_beta_lamda(phi_Nq = phi_Nqs.hi[2], 
                                  #f_Nq = f_Nqs[2], 
                                  #f_Nq = fNqs.ch[2],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3]
R0_fert_lo.hi = get_Ro_beta_lamda(phi_Nq = phi_Nqs.lo[2], 
                                  #f_Nq = f_Nqs[2],
                                  #f_Nq = fNqs.ch[2],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3]  

#atrazine+chlorpyrifos R0
R0_atch0.hi = get_Ro_beta_lamda(muPq = muPq, 
                                phi_Nq = phi_Nqs[3], 
                                #f_Nq = f_Nqs[3], 
                                #f_Nq = fNqs.ch[3],
                                beta = pt4$beta, 
                                lamda = pt4$lamda)[3]
R0_atch_up.hi = get_Ro_beta_lamda(muPq = muPq, 
                                  phi_Nq = phi_Nqs.hi[3], 
                                  #f_Nq = f_Nqs[3], 
                                  #f_Nq = fNqs.ch[3],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3]
R0_atch_lo.hi = get_Ro_beta_lamda(muPq = muPq, 
                                  phi_Nq = phi_Nqs.lo[3], 
                                  #f_Nq = f_Nqs[3], 
                                  #f_Nq = fNqs.ch[3],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3] 

#atrazine+fertilizer R0
R0_atfe0.hi = get_Ro_beta_lamda(phi_Nq = phi_Nqs[4], 
                                #f_Nq = f_Nqs[4], 
                                #f_Nq = fNqs.ch[4],
                                beta = pt4$beta, 
                                lamda = pt4$lamda)[3]
R0_atfe_up.hi = get_Ro_beta_lamda(phi_Nq = phi_Nqs.hi[4], 
                                  #f_Nq = f_Nqs[4], 
                                  #f_Nq = fNqs.ch[4],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3]
R0_atfe_lo.hi = get_Ro_beta_lamda(phi_Nq = phi_Nqs.lo[4], 
                                  #f_Nq = f_Nqs[4], 
                                  #f_Nq = fNqs.ch[4],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3] 

#chlorpyrifos+fertilizer R0
R0_chfe0.hi = get_Ro_beta_lamda(muPq = muPq, 
                                phi_Nq = phi_Nqs[2], 
                                #f_Nq = f_Nqs[2],
                                #f_Nq = fNqs.ch[2],
                                beta = pt4$beta, 
                                lamda = pt4$lamda)[3]
R0_chfe_up.hi = get_Ro_beta_lamda(muPq = muPq, 
                                  phi_Nq = phi_Nqs.hi[2], 
                                  #f_Nq = f_Nqs[2], 
                                  #f_Nq = fNqs.ch[2],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3]
R0_chfe_lo.hi = get_Ro_beta_lamda(muPq = muPq, 
                                  phi_Nq = phi_Nqs.lo[2], 
                                  #f_Nq = f_Nqs[2], 
                                  #f_Nq = fNqs.ch[2],
                                  beta = pt4$beta, 
                                  lamda = pt4$lamda)[3]  

#atrazine+chlorpyrifos+fertilizer R0
R0_tre0.hi = get_Ro_beta_lamda(muPq = muPq, 
                               phi_Nq = phi_Nqs[4], 
                               #f_Nq = f_Nqs[4], 
                               #f_Nq = fNqs.ch[4],
                               beta = pt4$beta, 
                               lamda = pt4$lamda)[3]
R0_tre_up.hi = get_Ro_beta_lamda(muPq = muPq, 
                                 phi_Nq = phi_Nqs.hi[4], 
                                 #f_Nq = f_Nqs[4], 
                                 #f_Nq = fNqs.ch[4],
                                 beta = pt4$beta, 
                                 lamda = pt4$lamda)[3]
R0_tre_lo.hi = get_Ro_beta_lamda(muPq = muPq, 
                                 phi_Nq = phi_Nqs.lo[4], 
                                 #f_Nq = f_Nqs[4], 
                                 #f_Nq = fNqs.ch[4],
                                 beta = pt4$beta, 
                                 lamda = pt4$lamda)[3] 

#Step 2.2 Plot R0 estimates across agrochemical treatment groups ##################
r0s.3<-data.frame("Treatment" = rep(c('At', 'Ch', 'Fe', 'At:Ch', 'At:Fe', 'Ch:Fe', 'At:Ch:Fe'), times=3),
                  "Season" = rep(c('Low', 'High', 'Mid'), each = 7),
                  'r0_0'= c(R0_atra0.low, R0_chlor0.low, R0_fert0.low, R0_atch0.low, R0_atfe0.low, R0_chfe0.low, R0_tre0.low,
                            R0_atra0.mid, R0_chlor0.mid, R0_fert0.mid, R0_atch0.mid, R0_atfe0.mid, R0_chfe0.mid, R0_tre0.mid,
                            R0_atra0.hi, R0_chlor0.hi, R0_fert0.hi, R0_atch0.hi, R0_atfe0.hi, R0_chfe0.hi, R0_tre0.hi),
                  'r0_up'= c(R0_atra_up.low, R0_chlor_up.low, R0_fert_up.low, R0_atch_up.low, R0_atfe_up.low, R0_chfe_up.low, R0_tre_up.low,
                             R0_atra_up.mid, R0_chlor_up.mid, R0_fert_up.mid, R0_atch_up.mid, R0_atfe_up.mid, R0_chfe_up.mid, R0_tre_up.mid,
                             R0_atra_up.hi, R0_chlor_up.hi, R0_fert_up.hi, R0_atch_up.hi, R0_atfe_up.hi, R0_chfe_up.hi, R0_tre_up.hi),
                  'r0_lo'= c(R0_atra_lo.low,  R0_chlor_lo.low, R0_fert_lo.low, R0_atch_lo.low, R0_atfe_lo.low, R0_chfe_lo.low, R0_tre_lo.low,
                             R0_atra_lo.mid,  R0_chlor_lo.mid, R0_fert_lo.mid, R0_atch_lo.mid, R0_atfe_lo.mid, R0_chfe_lo.mid, R0_tre_lo.mid,
                             R0_atra_lo.hi,  R0_chlor_lo.hi, R0_fert_lo.hi, R0_atch_lo.hi, R0_atfe_lo.hi, R0_chfe_lo.hi, R0_tre_lo.hi))  

r0s.3$Treatment<-factor(r0s.3$Treatment, levels = c('Fe', 'At', 'At:Fe', #Bottom up effects only
                                                    'Ch', #Top-down effects only
                                                    'Ch:Fe', 'At:Ch',   'At:Ch:Fe')) #Both top-down and bottom-up effects

r0s.3$Season<-factor(r0s.3$Season, levels = c('High', 'Mid', 'Low')) 

gg1<-ggplot(r0s.3, aes(x=Treatment, y=r0_0, group = Treatment, col = Season))+
  #Theme formatting
  theme_bw()+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=10))+
  scale_color_manual(values = c('red', 'black', 'blue')) +
  scale_y_continuous(breaks=c(0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0), limits=c(0,8))+
  xlab("")+
  ylab(expression('R'[0]))+
  #Village R0 lines
  geom_hline(aes(yintercept=R0_vil), colour='grey30', size=1, linetype=3)+
  #geom_hline(aes(yintercept=R0_up), colour='red', linetype=2, alpha=0.5, size=1)+
  #geom_hline(aes(yintercept=R0_lo), colour='red', linetype=2, alpha=0.5, size=1)+
  #geom_hline(aes(yintercept=R0_vil.p), colour='grey30', size=1, linetype=3)+
  #geom_hline(aes(yintercept=R0_up.p), colour='green2', linetype=2, alpha=0.5, size=1)+
  #geom_hline(aes(yintercept=R0_lo.p), colour='green2', linetype=2, alpha=0.5, size=1)+
  #Adding data from data frame
  geom_point(size=3) +
  geom_errorbar(aes(ymin=r0_lo,
                    ymax=r0_up),
                width=.1, position=position_dodge(.7))+
  #Add labels to R0 lines
  #geom_label(x=1.475, y=3.2, label="R ",  size=5, colour = 'grey30', fill='white', label.size = NA)+
  #geom_text(x=1.525, y=3.09, label="0,f",  size=3, colour = 'grey30')+
  #geom_label(x=5.475, y=0.53, label='R ',  size=5, colour = 'grey30', fill='white', label.size = NA)+
  #geom_text(x=5.525, y=0.508, label="0,f",  size=3, colour = 'grey30')+
  #geom_text(x=5.585, y=0.50, label="p",  size=3, colour = 'grey30')+
  #Add treatment labels
  geom_segment(x=0.7, xend=3.3, y=7, yend=7, colour='grey40', lineend='square')+
  geom_text(x=2, y=7.45, label='bottom-up effects', size=5, colour='grey40')+
  geom_segment(x=3.7, xend=4.3, y=7, yend=7, colour='grey40', lineend='square')+
  geom_text(x=4, y=7.65, label='top-down', size=5, colour='grey40')+
  geom_text(x=4, y=7.35, label='effects', size=5, colour='grey40')+
  geom_segment(x=4.7, xend=7.3, y=7, yend=7, colour='grey40', lineend='square')+
  geom_text(x=6, y=7.45, label='bottom-up & top-down effects', size=5, colour='grey40')+
  #Add plot label
  geom_text(label='A', x=0.57, y=8, size=10)

#Transmission parameter distributions##################
ftng<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_beta_2lamda_rand3.csv')

#Time-weighted average of lamda
ftng$lamda.twa<-(301/584)*ftng$lamda1+(253/584)*ftng$lamda2

#R0 using time weighted average
for(i in 1:nrow(ftng)){
  ftng[i,17]<-get_Ro_beta_lamda(muPq = p.dead, beta = ftng[i,1], lamda = ftng[i,16])[3]
}
colnames(ftng)[17]='R0.twa'

#Keep parameter combinations that fit model while maintaining reasonable snail prevalence
good<-subset(ftng, likelihood_range234 >0 & snail.prev <10)

plot(x=good$beta1, y=good$lamda1, pch = 16, xlab = 'beta', ylab = 'lamda', 
     ylim = c(0, max(good$lamda2)))
points(x = good$beta1, y = good$lamda2, pch = 16, col = 'red')
points(x = good$beta1, y = good$lamda.twa, pch = 16, col = 'blue')
legend('topright', legend = c('lamda1', 'lamda2', 'lamda.twa'), 
       pch = 16, col = c(1,2,4), cex=0.5)

hist(good$R0.twa) 

hist(good$beta1)
  mean.beta<-mean(good$beta1)
  sd.beta<-sd(good$beta1)
    hist(rnorm(10000, mean = mean.beta, sd = sd.beta))
hist(good$lamda1)
hist(good$lamda2)

hist(-log(good$likelihood_range234))

beta.best<-good$beta1[-log(good$likelihood_range234)==min(-log(good$likelihood_range234))]
  hist(rnorm(10000, mean = beta.best, sd = 0.2*beta.best))

lamda.best<-good$lamda.twa[-log(good$likelihood_range234)==min(-log(good$likelihood_range234))]
  hist(rnorm(10000, mean = lamda.best, sd = 0.2*lamda.best))


get_Ro_beta_lamda(muPq = p.dead, beta = beta.best, lamda = lamda.best)

#Generate probability of each parameter combination based on relative likelihood
for(i in 1:nrow(good)){
  good[i,18]=(-log(good[i,14])/sum(-log(good[,14])))
}

hist(good$V18)

#Sample distributions of parameters to use in R0 expression with probability ~ relative likelihood
beta.samp<-sample(good$beta1, size=1000000, replace = T, prob = good$V18)
  hist(beta.samp)

lamda.twa.samp<-sample(good$lamda.twa, size = 1000000, replace = T, prob = good$V18)
  hist(lamda.twa.samp)
  
#Updated transmission parameter distributions #################
  fin<-read.csv('fit.fin2.csv')
  
  fin<-subset(fin, snail.prev < 10 & R0 >1 & likelihood_range234 != 0)
  
  fin$likelihood_range234<- -log(fin$likelihood_range234)
  
  best<-unique(fin[which(fin$likelihood_range234 == min(fin$likelihood_range234)),])
  
  plot(x=fin$lamda.twa[fin$beta == best$beta], y = fin$likelihood_range234[fin$beta == best$beta])
  
    get_Ro_beta_lamda(muPq = 0, beta = best$beta, lamda = best$lamda.twa)
  
  samps<-subset(fin, select = c(beta, lamda1, lamda2, lamda.twa, likelihood_range234, prob)) 
    samps<-unique(samps) #Remove duplicates
    
    samps$prob = samps$likelihood_range234 / sum(samps$likelihood_range234)
    
  betasamp<-sample(samps$beta, size = 10000, replace = T, samps$prob)
    hist(betasamp)  
    
  ltwasamp<-sample(samps$lamda.twa, size = 10000, replace = T, samps$prob)
    hist(ltwasamp) 
    
  for(i in 1:length(betasamp)){
   R0samp[i] = get_Ro_beta_lamda(muPq = p.dead, phi_Nq = 1, beta = betasamp[i], lamda = ltwasamp[i])[3]
  }  
    hist(R0samp)
    
#Bottom-up parameter distributions ######################
  dat<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/R_use.csv")

#Bootstrapping to obtain distribution of Phi_Nq estimates in each bottom-up treatment 
    Fes<-dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==1]
    Ats<-dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==0]
    Fe.Ats<-dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==1]
    refs<-dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==0]
    ref.mean<-mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==0])
  
  
  Fe.samples <- lapply(1:100000, function(i)
    (sample(Fes, replace = T)/ref.mean))
      Fe.mean<-sapply(Fe.samples, mean)
      hist(Fe.mean, breaks = 15, xlim = c(0,3.5))
  
  At.samples <- lapply(1:100000, function(i)
    (sample(Ats, replace = T)/ref.mean))
      At.mean<-sapply(At.samples, mean)
      hist(At.mean, breaks = 15, xlim = c(0,3.5))
  
  FeAt.samples <- lapply(1:100000, function(i)
    (sample(Fe.Ats, replace = T)/ref.mean))
      FeAt.mean<-sapply(FeAt.samples, mean)
      hist(FeAt.mean, breaks = 15, xlim = c(0,3.5))
      
#Top-down parameter distributions ########################
  p.24<-dat$p.all_24[dat$chlor==0] #Prawn survivial in ChlorP free tanks
    p.24<-3-p.24 #convert to deaths
    #Bootstrap
    p.24.samples <- lapply(1:100000, function(i)
      (sample(p.24, replace = T)/3))
        muP.mean<-sapply(p.24.samples, mean)
        hist(muP.mean)
          abline(v=p.dead,col='red', lty=2)
        
  p.24.ch<-dat$p.all_24[dat$chlor==1] #Prawn survival in ChlorP pos tanks
    p.24.ch<-3-p.24.ch #convert to deaths
    #Bootstrap
    p.24.ch.samples <- lapply(1:100000, function(i)
      (sample(p.24.ch, replace = T)/3))
        muPq.mean<-sapply(p.24.ch.samples, mean)
        hist(muPq.mean)
          abline(v=p.dead,col='red', lty=2)
        
  
#R0 monte carlo simulations with samples from transmission and bottom-up parameters #######################
  parameters["mu_P"]<-0 #make all mortality correspond to observed mortalities in mesocosm
    
  control<-sapply(1:10000, 
                  function(i) get_Ro_beta_lamda(muPq = sample(muP.mean, 1),
                                                phi_Nq = 1,
                                                beta = sample(good$beta1, size=1, prob = good$V18),
                                                lamda = sample(good$lamda.twa, size=1, prob = good$V18))[3])       
    hist(control)
    mean(control)
    sd(control)       
    
  At.MC<-sapply(1:10000, 
                function(i) get_Ro_beta_lamda(muPq = sample(muP.mean, 1),
                                              phi_Nq = sample(At.mean, 1),
                                              beta = sample(good$beta1, size=1, prob = good$V18),
                                              lamda = sample(good$lamda.twa, size=1, prob = good$V18))[3])    
    hist(At.MC)
    mean(At.MC)
    sd(At.MC)
    
  
  Fe.MC<-sapply(1:10000, 
                function(i) get_Ro_beta_lamda(muPq = sample(muP.mean, 1),
                                              phi_Nq = sample(Fe.mean, 1),
                                              beta = sample(good$beta1, size=1, prob = good$V18),
                                              lamda = sample(good$lamda.twa, size=1, prob = good$V18))[3])    
    hist(Fe.MC)
    mean(Fe.MC)
    sd(Fe.MC)
    
  FeAt.MC<-sapply(1:10000, 
                  function(i) get_Ro_beta_lamda(muPq = sample(muP.mean, 1),
                                                phi_Nq = sample(FeAt.mean, 1),
                                                beta = sample(good$beta1, size=1, prob = good$V18),
                                                lamda = sample(good$lamda.twa, size=1, prob = good$V18))[3])    
    hist(FeAt.MC)
    mean(FeAt.MC)
    sd(FeAt.MC)
  
  Ch.MC<-sapply(1:10000, 
                function(i) get_Ro_beta_lamda(muPq = sample(muPq.mean, 1),
                                              phi_Nq = 1,
                                              beta = sample(good$beta1, size=1, prob = good$V18),
                                              lamda = sample(good$lamda.twa, size=1, prob = good$V18))[3])    
    hist(Ch.MC, xlim=c(0,15), ylim=c(0,2500))
    mean(Ch.MC)
    sd(Ch.MC)
    
  At.Ch.MC<-sapply(1:10000, 
                   function(i) get_Ro_beta_lamda(muPq = sample(muPq.mean, 1),
                                                 phi_Nq = sample(At.mean, 1),
                                                 beta = sample(good$beta1, size=1, prob = good$V18),
                                                 lamda = sample(good$lamda.twa, size=1, prob = good$V18))[3])    
    hist(At.Ch.MC, xlim=c(0,15), ylim=c(0,2500))
    mean(At.Ch.MC)
    sd(At.Ch.MC)
  
  Fe.Ch.MC<-sapply(1:10000, 
                   function(i) get_Ro_beta_lamda(muPq = sample(muPq.mean, 1),
                                                 phi_Nq = sample(Fe.mean, 1),
                                                 beta = sample(good$beta1, size=1, prob = good$V18),
                                                 lamda = sample(good$lamda.twa, size=1, prob = good$V18))[3])    
    hist(Fe.Ch.MC, xlim=c(0,15), ylim=c(0,2500))
    mean(Fe.Ch.MC)
    sd(Fe.Ch.MC)
  
  Tre.MC<-sapply(1:10000, 
                 function(i) get_Ro_beta_lamda(muPq = sample(muPq.mean, 1), 
                                               phi_Nq = sample(FeAt.mean, 1),
                                               beta = sample(good$beta1, size=1, prob = good$V18),
                                               lamda = sample(good$lamda.twa, size=1, prob = good$V18))[3])    
    hist(Tre.MC, xlim=c(0,15), ylim=c(0,2500))
    mean(Tre.MC)
    sd(Tre.MC)
#Violin plot ##################  
MC.df<-data.frame('R0' = c(At.MC, Fe.MC, FeAt.MC, 
                           Ch.MC, 
                           At.Ch.MC, Fe.Ch.MC, Tre.MC),
                  'treat' = c(rep('At', times=length(At.MC)), 
                              rep('Fe', times=length(Fe.MC)), 
                              rep('At:Fe', times=length(FeAt.MC)), 
                              rep('Ch', times=length(Ch.MC)), 
                              rep('At:Ch', times=length(At.Ch.MC)), 
                              rep('Ch:Fe', times=length(Fe.Ch.MC)), 
                              rep('At:Ch:Fe', times=length(Tre.MC))))

MC.df$treat<-factor(MC.df$treat, levels = c('Fe', 'At', 'At:Fe', #Bottom up effects only
                                            'Ch', #Top-down effects only
                                            'Ch:Fe', 'At:Ch', 'At:Ch:Fe')) #Both top-down and bottom-up effects

MC.df<-MC.df[complete.cases(MC.df),]

ggplot(MC.df, aes(x=treat, y=R0, group = treat))+
  #Theme formatting
  theme_bw()+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=10))+
  scale_y_continuous(breaks=seq(0, 15, 1), limits=c(0,max(MC.df$R0 + 1)))+
  xlab("")+
  ylab(expression('R'[0]))+
  #Village R0 lines
  #geom_hline(aes(yintercept=R0_vil), colour='grey30', size=1, linetype=3)+
  #geom_hline(aes(yintercept=R0_up), colour='red', linetype=2, alpha=0.5, size=1)+
  #geom_hline(aes(yintercept=R0_lo), colour='red', linetype=2, alpha=0.5, size=1)+
  #geom_hline(aes(yintercept=R0_vil.p), colour='grey30', size=1, linetype=3)+
  #geom_hline(aes(yintercept=R0_up.p), colour='green2', linetype=2, alpha=0.5, size=1)+
  #geom_hline(aes(yintercept=R0_lo.p), colour='green2', linetype=2, alpha=0.5, size=1)+
  #Adding data from data frame
  geom_violin()+
  stat_summary(fun.y=mean, geom='point', size=2)
  #geom_label(x=1.475, y=3.2, label="R ",  size=5, colour = 'grey30', fill='white', label.size = NA)+
  #geom_text(x=1.525, y=3.09, label="0,f",  size=3, colour = 'grey30')+
  #geom_label(x=5.475, y=0.53, label='R ',  size=5, colour = 'grey30', fill='white', label.size = NA)+
  #geom_text(x=5.525, y=0.508, label="0,f",  size=3, colour = 'grey30')+
  #geom_text(x=5.585, y=0.50, label="p",  size=3, colour = 'grey30')+
  #Add treatment labels
  geom_segment(x=0.7, xend=3.3, y=7, yend=7, colour='grey40', lineend='square')+
  geom_text(x=2, y=7.45, label='bottom-up effects', size=5, colour='grey40')+
  geom_segment(x=3.7, xend=4.3, y=7, yend=7, colour='grey40', lineend='square')+
  geom_text(x=4, y=7.65, label='top-down', size=5, colour='grey40')+
  geom_text(x=4, y=7.45, label='effects', size=5, colour='grey40')+
  geom_segment(x=4.7, xend=7.3, y=7, yend=7, colour='grey40', lineend='square')+
  geom_text(x=6, y=7.45, label='bottom-up & top-down effects', size=5, colour='grey40')+
  #Add plot label
  geom_text(label='A', x=0.57, y=8, size=10)
  
#Point/bar plot ##################
  MC.df2<-data.frame('treat' = c('control', 
                                 'Fe', 'At', 'At:Fe', #Bottom up effects only
                                 'Ch', #Top-down effects only
                                 'Ch:Fe', 'At:Ch', 'At:Ch:Fe'),
                     'R0.mean' = c(mean(control),
                                   mean(Fe.MC),
                                   mean(At.MC),
                                   mean(FeAt.MC),
                                   mean(Ch.MC),
                                   mean(Fe.Ch.MC),
                                   mean(At.Ch.MC),
                                   mean(Tre.MC)),
                     'R0.up' = c(mean(control) + 1.96*sd(control),
                                 mean(Fe.MC) + 1.96*sd(Fe.MC),
                                 mean(At.MC) + 1.96*sd(At.MC),
                                 mean(FeAt.MC) + 1.96*sd(FeAt.MC),
                                 mean(Ch.MC) + 1.96*sd(Ch.MC),
                                 mean(Fe.Ch.MC) + 1.96*sd(Fe.Ch.MC),
                                 mean(At.Ch.MC) + 1.96*sd(At.Ch.MC),
                                 mean(Tre.MC) + 1.96*sd(Tre.MC)),
                     'R0.lo' = c(mean(control) - 1.96*sd(control),
                                 mean(Fe.MC) - 1.96*sd(Fe.MC),
                                 mean(At.MC) - 1.96*sd(At.MC),
                                 mean(FeAt.MC) - 1.96*sd(FeAt.MC),
                                 mean(Ch.MC) - 1.96*sd(Ch.MC),
                                 mean(Fe.Ch.MC) - 1.96*sd(Fe.Ch.MC),
                                 mean(At.Ch.MC) - 1.96*sd(At.Ch.MC),
                                 mean(Tre.MC) - 1.96*sd(Tre.MC)))
  
  MC.df2$treat<-factor(MC.df2$treat, levels = c('control', 'Fe', 'At', 'At:Fe', #Bottom up effects only
                                              'Ch', #Top-down effects only
                                              'Ch:Fe', 'At:Ch', 'At:Ch:Fe')) #Both top-down and bottom-up effects
  MC.df2$R0.lo[MC.df2$R0.lo < 0]<-0
  
  ggplot(MC.df2, aes(x=treat, y=R0.mean, group = treat))+
    #Theme formatting
    theme_bw()+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=10))+
    scale_y_continuous(breaks=c(0,1,5,10), limits=c(0,max(MC.df2$R0.up + 3)))+
    xlab("")+
    ylab(expression('R'[0]))+
    #Village R0 lines
    #geom_hline(aes(yintercept=R0_vil), colour='grey30', size=1, linetype=3)+
    #geom_hline(aes(yintercept=R0_up), colour='red', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_lo), colour='red', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_vil.p), colour='grey30', size=1, linetype=3)+
    #geom_hline(aes(yintercept=R0_up.p), colour='green2', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_lo.p), colour='green2', linetype=2, alpha=0.5, size=1)+
    #Adding data from data frame
    geom_point(size=3) +
    geom_errorbar(aes(ymin=R0.lo,
                      ymax=R0.up),
                  width=.1, position=position_dodge(.7))+
  #geom_label(x=1.475, y=3.2, label="R ",  size=5, colour = 'grey30', fill='white', label.size = NA)+
  #geom_text(x=1.525, y=3.09, label="0,f",  size=3, colour = 'grey30')+
  #geom_label(x=5.475, y=0.53, label='R ',  size=5, colour = 'grey30', fill='white', label.size = NA)+
  #geom_text(x=5.525, y=0.508, label="0,f",  size=3, colour = 'grey30')+
  #geom_text(x=5.585, y=0.50, label="p",  size=3, colour = 'grey30')+
  #Add treatment labels
  geom_segment(x=1.7, xend=4.3, y=11, yend=11, colour='grey40', lineend='square')+
    geom_text(x=3, y=11.45, label='bottom-up effects', size=5, colour='grey40')+
    geom_segment(x=4.7, xend=5.3, y=11, yend=11, colour='grey40', lineend='square')+
    geom_text(x=5, y=11.95, label='top-down', size=5, colour='grey40')+
    geom_text(x=5, y=11.45, label='effects', size=5, colour='grey40')+
    geom_segment(x=5.7, xend=8.3, y=11, yend=11, colour='grey40', lineend='square')+
    geom_text(x=7, y=11.45, label='bottom-up & top-down effects', size=5, colour='grey40')+
    #Add plot label
    geom_text(label='A', x=0.57, y=8, size=10)
  

#R0 monte carlo simulations with mesocosm uncertainty ONLY #######################
  parameters["mu_P"]<-0 #make all mortality correspond to observed mortalities in mesocosm
  
  control<-sapply(1:10000, 
                  function(i) get_Ro_beta_lamda(muPq = sample(muP.mean, 1),
                                                phi_Nq = 1,
                                                beta = beta.best,
                                                lamda = lamda.best)[3])       
  hist(control)
  mean(control)
  sd(control)       
  
  At.MC<-sapply(1:10000, 
                function(i) get_Ro_beta_lamda(muPq = sample(muP.mean, 1),
                                              phi_Nq = sample(At.mean, 1),
                                              beta = beta.best,
                                              lamda = lamda.best)[3])    
  hist(At.MC)
  mean(At.MC)
  sd(At.MC)
  
  
  Fe.MC<-sapply(1:10000, 
                function(i) get_Ro_beta_lamda(muPq = sample(muP.mean, 1),
                                              phi_Nq = sample(Fe.mean, 1),
                                              beta = beta.best,
                                              lamda = lamda.best)[3])    
  hist(Fe.MC)
  mean(Fe.MC)
  sd(Fe.MC)
  
  FeAt.MC<-sapply(1:10000, 
                  function(i) get_Ro_beta_lamda(muPq = sample(muP.mean, 1),
                                                phi_Nq = sample(FeAt.mean, 1),
                                                beta = beta.best,
                                                lamda = lamda.best)[3])    
  hist(FeAt.MC)
  mean(FeAt.MC)
  sd(FeAt.MC)
  
  Ch.MC<-sapply(1:10000, 
                function(i) get_Ro_beta_lamda(muPq = sample(muPq.mean, 1),
                                              phi_Nq = 1,
                                              beta = beta.best,
                                              lamda = lamda.best)[3])    
  hist(Ch.MC)
  mean(Ch.MC)
  sd(Ch.MC)
  
  At.Ch.MC<-sapply(1:10000, 
                   function(i) get_Ro_beta_lamda(muPq = sample(muPq.mean, 1),
                                                 phi_Nq = sample(At.mean, 1),
                                                 beta = beta.best,
                                                 lamda = lamda.best)[3])    
  hist(At.Ch.MC)
  mean(At.Ch.MC)
  sd(At.Ch.MC)
  
  Fe.Ch.MC<-sapply(1:10000, 
                   function(i) get_Ro_beta_lamda(muPq = sample(muPq.mean, 1),
                                                 phi_Nq = sample(Fe.mean, 1),
                                                 beta = beta.best,
                                                 lamda = lamda.best)[3])    
  hist(Fe.Ch.MC)
  mean(Fe.Ch.MC)
  sd(Fe.Ch.MC)
  
  Tre.MC<-sapply(1:10000, 
                 function(i) get_Ro_beta_lamda(muPq = sample(muPq.mean, 1), 
                                               phi_Nq = sample(FeAt.mean, 1),
                                               beta = beta.best,
                                               lamda = lamda.best)[3])    
  hist(Tre.MC)
  mean(Tre.MC)
  sd(Tre.MC)
  
MC.df3<-data.frame('R0' = c(At.MC, Fe.MC, FeAt.MC, 
                             Ch.MC, 
                             At.Ch.MC, Fe.Ch.MC, Tre.MC),
                    'treat' = c(rep('At', times=length(At.MC)), 
                                rep('Fe', times=length(Fe.MC)), 
                                rep('At:Fe', times=length(FeAt.MC)), 
                                rep('Ch', times=length(Ch.MC)), 
                                rep('At:Ch', times=length(At.Ch.MC)), 
                                rep('Ch:Fe', times=length(Fe.Ch.MC)), 
                                rep('At:Ch:Fe', times=length(Tre.MC))))
  
MC.df3$treat<-factor(MC.df$treat, levels = c('Fe', 'At', 'At:Fe', #Bottom up effects only
                                              'Ch', #Top-down effects only
                                              'Ch:Fe', 'At:Ch', 'At:Ch:Fe')) #Both top-down and bottom-up effects
  
  ggplot(MC.df3, aes(x=treat, y=R0, group = treat))+
    #Theme formatting
    theme_bw()+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=10))+
    scale_y_continuous(breaks=seq(0, 15, 1), limits=c(0,max(MC.df$R0 + 3)))+
    xlab("")+
    ylab(expression('R'[0]))+
    #Village R0 lines
    #geom_hline(aes(yintercept=R0_vil), colour='grey30', size=1, linetype=3)+
    #geom_hline(aes(yintercept=R0_up), colour='red', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_lo), colour='red', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_vil.p), colour='grey30', size=1, linetype=3)+
    #geom_hline(aes(yintercept=R0_up.p), colour='green2', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_lo.p), colour='green2', linetype=2, alpha=0.5, size=1)+
    #Adding data from data frame
    geom_violin(width=1)

#Point/bar plot ##################
  MC.df4<-data.frame('treat' = c('Fe', 'At', 'At:Fe', #Bottom up effects only
                                 'Ch', #Top-down effects only
                                 'Ch:Fe', 'At:Ch', 'At:Ch:Fe'),
                     'R0.mean' = c(mean(Fe.MC),
                                   mean(At.MC),
                                   mean(FeAt.MC),
                                   mean(Ch.MC),
                                   mean(Fe.Ch.MC),
                                   mean(At.Ch.MC),
                                   mean(Tre.MC)),
                     'R0.up' = c(mean(Fe.MC) + 1.96*sd(Fe.MC),
                                 mean(At.MC) + 1.96*sd(At.MC),
                                 mean(FeAt.MC) + 1.96*sd(FeAt.MC),
                                 mean(Ch.MC) + 1.96*sd(Ch.MC),
                                 mean(Fe.Ch.MC) + 1.96*sd(Fe.Ch.MC),
                                 mean(At.Ch.MC) + 1.96*sd(At.Ch.MC),
                                 mean(Tre.MC) + 1.96*sd(Tre.MC)),
                     'R0.lo' = c(mean(Fe.MC) - 1.96*sd(Fe.MC),
                                 mean(At.MC) - 1.96*sd(At.MC),
                                 mean(FeAt.MC) - 1.96*sd(FeAt.MC),
                                 mean(Ch.MC) - 1.96*sd(Ch.MC),
                                 mean(Fe.Ch.MC) - 1.96*sd(Fe.Ch.MC),
                                 mean(At.Ch.MC) - 1.96*sd(At.Ch.MC),
                                 mean(Tre.MC) - 1.96*sd(Tre.MC)))
  
  MC.df4$treat<-factor(MC.df4$treat, levels = c('Fe', 'At', 'At:Fe', #Bottom up effects only
                                                'Ch', #Top-down effects only
                                                'Ch:Fe', 'At:Ch', 'At:Ch:Fe')) #Both top-down and bottom-up effects
  MC.df4$R0.lo[MC.df2$R0.lo < 0]<-0
  
  ggplot(MC.df4, aes(x=treat, y=R0.mean, group = treat))+
    #Theme formatting
    theme_bw()+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=10))+
    scale_y_continuous(breaks=c(0,1,5,10), limits=c(0,max(MC.df2$R0.up + 3)))+
    xlab("")+
    ylab(expression('R'[0]))+
    #Village R0 lines
    #geom_hline(aes(yintercept=R0_vil), colour='grey30', size=1, linetype=3)+
    #geom_hline(aes(yintercept=R0_up), colour='red', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_lo), colour='red', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_vil.p), colour='grey30', size=1, linetype=3)+
    #geom_hline(aes(yintercept=R0_up.p), colour='green2', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_lo.p), colour='green2', linetype=2, alpha=0.5, size=1)+
    #Adding data from data frame
    geom_point(size=3) +
    geom_errorbar(aes(ymin=R0.lo,
                      ymax=R0.up),
                  width=.1, position=position_dodge(.7))+
    #geom_label(x=1.475, y=3.2, label="R ",  size=5, colour = 'grey30', fill='white', label.size = NA)+
    #geom_text(x=1.525, y=3.09, label="0,f",  size=3, colour = 'grey30')+
    #geom_label(x=5.475, y=0.53, label='R ',  size=5, colour = 'grey30', fill='white', label.size = NA)+
    #geom_text(x=5.525, y=0.508, label="0,f",  size=3, colour = 'grey30')+
    #geom_text(x=5.585, y=0.50, label="p",  size=3, colour = 'grey30')+
    #Add treatment labels
    geom_segment(x=0.7, xend=3.3, y=11, yend=11, colour='grey40', lineend='square')+
    geom_text(x=2, y=11.45, label='bottom-up effects', size=5, colour='grey40')+
    geom_segment(x=3.7, xend=4.3, y=11, yend=11, colour='grey40', lineend='square')+
    geom_text(x=4, y=11.95, label='top-down', size=5, colour='grey40')+
    geom_text(x=4, y=11.45, label='effects', size=5, colour='grey40')+
    geom_segment(x=4.7, xend=7.3, y=11, yend=11, colour='grey40', lineend='square')+
    geom_text(x=6, y=11.45, label='bottom-up & top-down effects', size=5, colour='grey40')+
    #Add plot label
    geom_text(label='A', x=0.57, y=8, size=10)
  
  