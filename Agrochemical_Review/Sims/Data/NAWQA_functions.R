require(tidyverse)

#Functions for NAWQA database ##############
#Function to return NAWQA data matching a particular chemical name
get_nawqa_dat <- function(chem_name){
  dat <- pestdat %>% 
    filter(grepl(chem_name, LONGNAME, ignore.case = T)) %>% 
    mutate(mug_perL = if_else(grepl("micrograms per liter", LONGNAME), 1, 0),
           CONCENTRATION = if_else(mug_perL == 1, CONCENTRATION, CONCENTRATION/1000)) %>% 
    select(CONSTIT, LONGNAME, CONCENTRATION, REMARK, SITE_TYPE, SAMPLE_DATE)
  
  return(dat)
}

#Function to return max value in NAWQA data for particular chemical
get_nawqa_max <- function(chem_name){
    dat <- pestdat %>% 
      filter(grepl(chem_name, LONGNAME, ignore.case = T)) %>% 
      mutate(mug_perL = if_else(grepl("micrograms per liter", LONGNAME), 1, 0),
             CONCENTRATION = if_else(mug_perL == 1, CONCENTRATION, CONCENTRATION/1000)) %>% 
      pull(CONCENTRATION)
  
  if(length(dat) == 0){
    return(NA)
  } else  {
    return(max(as.numeric(dat)))
  }
    
}

#Function to return median and IQR in NAWQA for each chem
get_nawqa_sum <- function(chem_name, all_obs = FALSE){
  dat <- pestdat %>% filter(grepl(chem_name, LONGNAME, ignore.case = T)) %>% 
    mutate(mug_perL = if_else(grepl("micrograms per liter", LONGNAME), 1, 0),
           CONCENTRATION = if_else(mug_perL == 1, CONCENTRATION, CONCENTRATION/1000))
    
    #All observations
    all_obs <- as.numeric(dat %>% pull(CONCENTRATION))
    #Observations above detectable limit
    obs_obs <- as.numeric(dat %>% filter(REMARK != "<") %>% pull(CONCENTRATION))
    
    all_sum <- c(median(all_obs), quantile(all_obs, 0.25), quantile(all_obs, 0.75))
    obs_sum <- c(median(obs_obs), quantile(obs_obs, 0.25), quantile(obs_obs, 0.75))
  
#Return summary of all observations or just those that were above detection threshold?    
  if(all_obs){
    return(obs_sum)
  } else {
    return(all_sum)
  }
  
}

#Functions for midwest 2013 database ##############
get_midwest_pest2013_dat <- function(chem_name){
  dat <- midwest_pest_2013_fin %>% 
    filter(grepl(chem_name, PARM_NM, ignore.case = T))
  
  return(dat)
}

#Function to return max value in NAWQA data for particular chemical
get_midwest_pest2013_max <- function(chem_name){
    dat <- midwest_pest_2013_fin %>% 
      filter(grepl(chem_name, PARM_NM, ignore.case = T)) %>% 
      pull(RESULT_CEN)
  
  if(length(dat) == 0){
    return(NA)
  } else  {
    return(max(as.numeric(dat)))
  }

}

#Function to return median and IQR in NAWQA for each chem
get_midwest_pest2013_sum <- function(chem_name, all_obs = FALSE){
    dat <- midwest_pest_2013_fin %>% filter(grepl(chem_name, PARM_NM, ignore.case = T)) 
    
    #All observations
    all_obs <- as.numeric(dat %>% pull(RESULT_CEN))
    #Observations above detectable limit
    obs_obs <- as.numeric(dat %>% filter(REMARK_CEN != "<") %>% pull(RESULT_CEN))
    
    all_sum <- c(median(all_obs), quantile(all_obs, 0.25), quantile(all_obs, 0.75))
    obs_sum <- c(median(obs_obs), quantile(obs_obs, 0.25), quantile(obs_obs, 0.75))
  
  if(all_obs){
    return(obs_sum)
  } else {
    return(all_sum)
  }
  
}

