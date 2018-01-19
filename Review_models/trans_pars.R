#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

fin<-read.csv("~/RemaisWork/Schisto/R Codes/ag_schist/Gen1/shortlist_transParams.csv")
fin100 = subset(fin, negLL > sort(fin$negLL, decreasing = T)[101])  

  fin$weight = fin$negLL / sum(fin$negLL)
  
  fin.lam = sum(fin$lamda.twa * fin$weight)
  fin.lam1 = sum(fin$lamda1 * fin$weight)
  fin.lam2 = sum(fin$lamda2 * fin$weight)
  
  fin.beta = sum(fin$beta * fin$weight)
  
fin100$weight = fin100$negLL / sum(fin100$negLL)

  fin100.lam = sum(fin100$lamda.twa * fin100$weight)
  fin100.lam1 = sum(fin100$lamda1 * fin100$weight)
  fin100.lam2 = sum(fin100$lamda2 * fin100$weight)
  
  fin100.beta = sum(fin100$beta * fin100$weight)