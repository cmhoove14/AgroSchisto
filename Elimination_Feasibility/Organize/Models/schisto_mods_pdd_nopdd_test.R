#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

##################################################################################################

require(deSolve)
require(graphics)
require(tidyverse)

source("Elimination_Feasibility/Organize/Models/Reff_BBR_fns.R")
source("Elimination_Feasibility/Organize/Models/schisto_mods_pdd_nopdd.R")

#Keep low transmission parameter sets
outputfile<-"Elimination_Feasibility/Organize/Models/Outputs_Refs/shortlist_51.csv"
  shortlist_first100<-as.data.frame(read.csv(file=outputfile, header=TRUE, sep=","))
  sort_ind<-order(shortlist_first100$R0, decreasing=FALSE)
  shortlist_first100<-shortlist_first100[sort_ind,]
    shortlist_first100 = shortlist_first100[c(1:100),] #Get rid of other parameter estimates
  
time<-seq(from=0, to=40*365, by=2)
nstart=c(S=5000,E=0,I=0, Wt=10, Wu=10)
  
  output<-as.data.frame(ode(nstart,time,schisto_mod_nopdd,pars_Chris1))
  output.addNDD<-as.data.frame(ode(nstart,time,schisto_mod_nopdd_add_ndds,pars_Chris1))
  output.PDD<-as.data.frame(ode(nstart,time,schisto_mod_pdd,pars_Chris1))
  output.PDD.addNDD<-as.data.frame(ode(nstart,time,schisto_mod_pdd_add_ndds,pars_Chris1))

plot_sim_worms <- function(output.df){
  output.df %>%
    select(time, Wt, Wu) %>%
    gather(key = "Treatment", value = "Worm_Burden", -time) %>%
    ggplot(aes(x = time, y = Worm_Burden, lty = Treatment)) + geom_line(col = "purple")
  
}

plot_sim_snails <- function(output.df){
  output.df %>%
    select(time, S, E, I) %>%
    gather(key = "Infection_Class", value = "Density", -time) %>%
    ggplot(aes(x = time, y = Density, col = Infection_Class)) + geom_line()

}

plot_sim_worms(output)
plot_sim_worms(output.addNDD)
plot_sim_worms(output.PDD)
plot_sim_worms(output.PDD.addNDD)

plot_sim_snails(output)
plot_sim_snails(output.addNDD)
plot_sim_snails(output.PDD)
plot_sim_snails(output.PDD.addNDD)
