#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/Hasheesh2011_chlorpyrifos_profenofos_snails_fit.R")

#Hasheesh and Mohamed 2011 data and analysis assessing toxicity of Chlorpyrifos and Profenofos ###############
# Direct toxicity to snails (Bu. truncatus); daily mortality rate ########

#Chlorpyrifos #########
#Plot
mun.hsh.ch = data.frame(conc = c(0,.72, 1.32, 2.82)*1000,
                        mort = c(0,.25, .50 , .90),
                        surv = 0)
  mun.hsh.ch$surv = 1 - mun.hsh.ch$mort

png("Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Snails/hasheesh2011_chlorpyrifos_snail_mortality.png")

  plot(mun.hsh.ch$conc, mun.hsh.ch$mort, pch = 16, ylim = c(-0.1,1), xlim = c(0,3500),
       xlab = 'Chlorpyrifos (ppb)', ylab = 'prop dead')
    segments(x0 = 880, y0 = 0.5, x1 = 1980, y1 = 0.5)
  
    set.seed(43093)
    
  points(seq(0,3500,10), sapply(seq(0,3500,10), muNq_ch_hash11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
 
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  
  

#Profenofos #########
#Plot  
mun.hsh.prof = data.frame(conc = c(0,1.4, 2.5, 3.72)*1000,
                    mort = c(0,.25, .50 , .90),
                    surv = 0)
  mun.hsh.prof$surv = 1 - mun.hsh.prof$mort

png("Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Snails/hasheesh2011_profenofos_snail_mortality.png")
    
  plot(mun.hsh.prof$conc, mun.hsh.prof$mort, pch = 16, ylim = c(0,1), xlim = c(0,5000),
       xlab = 'Profenofos (ppm)', ylab = 'prop dead')
    segments(x0 = 1880, y0 = 0.5, x1 = 3330, y1 = 0.5)
  
    set.seed(43093)
    
  points(seq(0,5000,10), sapply(seq(0,5000,10), muNq_prof_hash11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  

#Direct toxicity to snails affecting reproduction (Table 2) #######
#longitudinal measurement of mean eggs/snail          
png("Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Snails/hasheesh2011_longitudinal_snail_reproduction.png")

  plot(fn$time[fn$chem == 'control'], fn$eggs_snail_day[fn$chem == 'control'], type = 'l', lwd=2,
        xlab = 'time (weeks)', ylab = 'eggs/snail/day', ylim = c(0,max(fn$eggs_snail_day)+2),
        main = 'Snail reproduction over time Hasheesh 2011')
    lines(fn$time[fn$chem == 'chlorpyrifos'], fn$eggs_snail_day[fn$chem == 'chlorpyrifos'], col = 'red', lwd=2)
    lines(fn$time[fn$chem == 'profenofos'], fn$eggs_snail_day[fn$chem == 'profenofos'], col = 'orange', lwd=2) 
    legend('topleft', legend = c('control', 'chlorpyrifos', 'profenofos'),
           lwd = 2, col = c(1,2, 'orange'), cex = 0.8, bty = 'n')
dev.off()    

#Sampling uncertainty functions in control and single dose groups
set.seed(43093)
fn_hash_ctrl_samp <- rnorm(10000, 1, 0.25) #just for comparison
fn_hash_ch_samp <- replicate(10000, fNq.hash.chlor.uncertainty())
fn_hash_prof_samp <- replicate(10000, fNq.hash.prof.uncertainty())

png("Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Snails/hasheesh2011_snail_reproduction_compare.png")

plot(density(fn_hash_ctrl_samp), lwd = 2, main = "f_N estimate from longitudinal snail reproduction data",
     xlim = c(-0.5,2), ylim = c(0,5.5))
    lines(density(fn_hash_ch_samp), col = 'red', lwd=2)
    lines(density(fn_hash_prof_samp), col = 'orange', lwd=2) 
    legend('topright', legend = c('control', 'chlorpyrifos', 'profenofos'),
           lwd = 2, col = c(1,2, 'orange'), cex = 0.8, bty = 'n')
    
dev.off()    
