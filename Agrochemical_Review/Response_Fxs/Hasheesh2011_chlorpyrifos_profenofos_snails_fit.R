#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Models/litchfield_wilcoxon_get_b1_from_slope.R")

#Hasheesh and Mohamed 2011 data and analysis assessing toxicity of Chlorpyrifos and Profenofos ###############
# Direct toxicity to snails (Bu. truncatus); daily mortality rate ########

#Chlorpyrifos #########
  lc50.ch.hash.report = 1.32
  slp.ch.hash.report = 2.5
  b1.ch.hash.report = get_b1(slp.ch.hash.report)
    se.lc50.ch.hash = mean(c(log10(1.98/lc50.ch.hash.report), 
                             log10(lc50.ch.hash.report/0.88))) / 1.96
  
  muNq_ch_hash11_uncertainty = function(In){
    Ins = (In/1000)
    lc50 = 10^(rnorm(1, log10(lc50.ch.hash.report), se.lc50.ch.hash))
    mun = pnorm(b1.ch.hash.report * log10(Ins/lc50)) 

    return(mun)
  }
  
#keep vector
keep.hsh.ch = c('muNq_ch_hash11_uncertainty', 'lc50.ch.hash.report', 'se.lc50.ch.hash', 'b1.ch.hash.report')    
  
#Profenofos #########
  lc50.prof.hash.report = 2.5
  slp.prof.hash.report = 1.6
  b1.prof.hash.report = get_b1(slp.prof.hash.report)
    se.lc50.prof.hash = mean(c(log10(3.33/lc50.prof.hash.report), log10(lc50.prof.hash.report/1.88))) / 1.96
  
  muNq_prof_hash11_uncertainty = function(In){
    Ins = (In/1000)
    lc50 = 10^(rnorm(1, log10(lc50.prof.hash.report), se.lc50.prof.hash))
    mun = pnorm((b1.prof.hash.report) * log10(Ins/lc50))

    return(mun)
  }
  
#keep vector
keep.hsh.prof = c('muNq_prof_hash11_uncertainty', 'lc50.prof.hash.report', 'se.lc50.prof.hash', 'b1.prof.hash.report')    

#Direct toxicity to snails affecting reproduction (Table 2) #######
fn<-read.csv('Agrochemical_Review/Response_Fxs/Data/Hasheesh2011_snail_mort_repro_weekly.csv')

#get estimate of mean eggs/snail/day for each treatment
  ctrl.eggs = fn$eggs_snail_day[fn$chem == 'control' & fn$surv != 0][-1]

  chlor.eggs = fn$eggs_snail_day[fn$chem == 'chlorpyrifos' & fn$surv != 0][-1]

  prof.eggs = fn$eggs_snail_day[fn$chem == 'profenofos' & fn$surv != 0][-1]

#Functions that estimate parameter value in single concentration group #########
    fN.hash.ctrl.uncertainty = function(){ #function based on snail egg masses
      ctrl1 = sample(ctrl.eggs, 1)
      ctrl2 = sample(ctrl.eggs, 1)
      fN = ctrl2 / ctrl1
      
      return(fN)
    }
  
    fNq.hash.chlor.uncertainty = function(){ #function based on snail egg masses
      ctrl = sample(ctrl.eggs, 1)
      chlor = sample(chlor.eggs, 1)
      fN = chlor / ctrl
      
      return(fN)
    }
    
    fNq.hash.prof.uncertainty = function(){ #function based on snail hatchlings
      ctrl = sample(ctrl.eggs, 1)
      prof = sample(prof.eggs, 1)
      fN = prof / ctrl
      
      return(fN)
    }  
    
keep.hsh.ch = c(keep.hsh.ch, 'fNq.hash.chlor.uncertainty', 'chlor.eggs', 'ctrl.eggs')

keep.hsh.prof = c(keep.hsh.prof, 'fNq.hash.prof.uncertainty', 'prof.eggs', 'ctrl.eggs')

keep.hsh.all = c(keep.hsh.ch, keep.hsh.prof)  
