#Get pesticide monitoring data from USGS website <https://www.sciencebase.gov/catalog/item/595e6b16e4b0d1f9f05702fd>

require(data.table)

pestsamp <- fread("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__1f%2F1b%2Ff1%2F1f1bf16aadeab0625b4ecc967d01f1dba0b80300", skip = 12)

  pestsamp_meta <- read.csv(url("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__1f%2F1b%2Ff1%2F1f1bf16aadeab0625b4ecc967d01f1dba0b80300"), nrows = 12)

pestname <- fread("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__95%2F9c%2F6c%2F959c6c9b8343a3ea72ae2384a76a8359ff2ea08d", skip = 6)

  pestname_meta <- read.csv(url("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__95%2F9c%2F6c%2F959c6c9b8343a3ea72ae2384a76a8359ff2ea08d"), nrows = 6)

pestsites <- fread("https://www.sciencebase.gov/catalog/file/get/595e6b16e4b0d1f9f05702fd?f=__disk__a2%2F8e%2Fd2%2Fa28ed208752341daeaa642ae562657ed00df7a12", skip = 13)

pestdat <- pestsamp %>% 
  inner_join(pestname %>% select(-CONSTIT), by = "PARM_CD") %>% 
  inner_join(pestsites %>% select(-SITE_QW_ID), by = "SITE_ABB")

#Get response function summary to access relevant agrochemicals
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

chems <- unique(rfx_sum$Chemical)

#Function to return NAWQA data matching a particular chemical name
get_nawqa_dat <- function(chem_name){
  dat <- pestdat %>% 
    filter(grepl(chem_name, LONGNAME, ignore.case = T)) %>% 
    filter(REMARK != "<") %>% select(CONCENTRATION,CONSTIT, LONGNAME, SITE_TYPE)
  
  return(dat)
}