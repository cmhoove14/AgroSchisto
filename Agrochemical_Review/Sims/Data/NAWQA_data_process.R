#Get pesticide monitoring data from USGS website <https://www.sciencebase.gov/catalog/item/595e6b16e4b0d1f9f05702fd>

require(data.table)

pestsamp <- fread("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__1f%2F1b%2Ff1%2F1f1bf16aadeab0625b4ecc967d01f1dba0b80300", skip = 12)

  pestsamp_meta <- read.csv(url("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__1f%2F1b%2Ff1%2F1f1bf16aadeab0625b4ecc967d01f1dba0b80300"), nrows = 12)

pestname <- fread("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__95%2F9c%2F6c%2F959c6c9b8343a3ea72ae2384a76a8359ff2ea08d", skip = 6)

  pestname_meta <- read.csv(url("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__95%2F9c%2F6c%2F959c6c9b8343a3ea72ae2384a76a8359ff2ea08d"), nrows = 6)

pestsites <- fread("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__a2%2F8e%2Fd2%2Fa28ed208752341daeaa642ae562657ed00df7a12", skip = 13)

#Glyphosate data obtained from different source since not included in NAWQA master spreadsheet
gly_dat <- fread("https://www.sciencebase.gov//catalog/file/get/5928a5dce4b016f7a93f8d7b?f=__disk__f7%2Ff0%2F8e%2Ff7f08ec5b316028dbe466c78e1159a196d1101db", skip = 1)

pestdat <- pestsamp %>% 
  inner_join(pestname %>% select(-CONSTIT), by = "PARM_CD") %>% 
  inner_join(pestsites %>% select(-SITE_QW_ID), by = "SITE_ABB")

#Get response function summary to access relevant agrochemicals
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

chems <- unique(rfx_sum$Chemical)

#Function to return NAWQA data matching a particular chemical name
get_nawqa_dat <- function(chem_name){
  if(chem_name %in% c("glyphosate", "Glyphosate")){
    dat <- gly_dat %>% 
      select(SAMPLE_DATE, GLYPHOSATE_RMK_ELISA, GLYPHOSATE_ELISA_ng_L) %>% 
      mutate(GLYPHOSATE_ELISA_ug_L = GLYPHOSATE_ELISA_ng_L/1000)
    
  } else {
    
    dat <- pestdat %>% 
      filter(grepl(chem_name, LONGNAME, ignore.case = T)) %>% 
      select(CONCENTRATION, CONSTIT, LONGNAME, REMARK, SITE_TYPE)
  
  }
  
  return(dat)
}

#Function to return max value in NAWQA data for particular chemical
get_nawqa_max <- function(chem_name){
  if(chem_name %in% c("glyphosate", "Glyphosate")){
    dat <- gly_dat %>% 
      mutate(GLYPHOSATE_ELISA_ug_L = GLYPHOSATE_ELISA_ng_L/1000) %>% 
      pull(GLYPHOSATE_ELISA_ug_L)
    
  } else {
    
    dat <- pestdat %>% 
      filter(grepl(chem_name, LONGNAME, ignore.case = T)) %>% 
      pull(CONCENTRATION)
  
  }
  
  return(max(as.numeric(dat)))
}

#Function to return range to use in r0 simulations
get_range <- function(nawqa_vals, peak_eec){
  max_ac <- max(c(nawqa_vals, peak_eec))
  
  ac_range <- c(0, 
                exp(seq(log(min(nawqa_vals)),
                        log(max_ac),
                        length.out = 100)),
                seq(ceiling(max_ac+0.01),
                    ceiling(max_ac+0.01)*1.25,
                    length.out = 10))
  
  return(ac_range)
}

save.image("~/RemaisWork/Schisto/R Codes/ag_schist/Agrochemical_Review/Sims/Data/NAWQA_dat_functions.RData")