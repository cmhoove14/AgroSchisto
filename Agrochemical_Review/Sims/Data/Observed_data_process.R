#Get pesticide monitoring data from USGS website <https://www.sciencebase.gov/catalog/item/595e6b16e4b0d1f9f05702fd>

require(data.table)
require(lubridate)
require(tidyverse)

#NAWQA Data ###########
pestsamp <- fread("https://www.sciencebase.gov/catalog/file/get/5af49c2ae4b0da30c1b44e2b?f=__disk__94%2F66%2Fc7%2F9466c74b9a0f3b509950d893ee0b310455fec3ad", skip = 12)

pestname <- fread("https://www.sciencebase.gov/catalog/file/get/5af49c2ae4b0da30c1b44e2b?f=__disk__0f%2F03%2F30%2F0f0330574a9a25b1712677efe64ce9cc565dfe9f", skip = 6)

pestsites <- fread("https://www.sciencebase.gov/catalog/file/get/5af49c2ae4b0da30c1b44e2b?f=__disk__24%2F2e%2Fb2%2F242eb2ee19e8c2f6c70e0aafe7071638159ed146", skip = 13)

pestmetadata <- xmlParse(fread("https://www.sciencebase.gov/catalog/file/get/5af49c2ae4b0da30c1b44e2b?f=__disk__cf%2F5d%2Fe8%2Fcf5de87bd66ea59184465d87c0ae56dd8689ed7e"))

pestdat <- pestsamp %>% 
  inner_join(pestname %>% dplyr::select(-CONSTIT), by = "PARM_CD") %>% 
  inner_join(pestsites %>% dplyr::select(-SITE_QW_ID), by = "SITE_ABB") %>% 
  mutate(SAMPLE_DATE = ymd(str_split(DATETIME, " ", simplify = TRUE)[,1]))

saveRDS(pestdat, "Agrochemical_Review/Sims/Data/NAWQA_pesticides.rds")

#Another USGS source for pesticide data from monitoring of midwest Ag streams in 2013: The NAWQA Midwest Stream Quality Assessment (MSQA) <https://www.sciencebase.gov/catalog/item/5928a5dce4b016f7a93f8d7b> ###############

msqa2013_dictionary <- fread("https://www.sciencebase.gov/catalog/file/get/5928a5dce4b016f7a93f8d7b?f=__disk__8f%2F4b%2F6e%2F8f4b6e114731e7efa99e2a1455d8b46e2ac47ef9", skip = 1)

msqa2013_all_dat <- fread("https://www.sciencebase.gov/catalog/file/get/5928a5dce4b016f7a93f8d7b?f=__disk__81%2Fab%2Fe2%2F81abe26219aba700c4d9bce0419812ff439128f4", skip = 1)

msqa2013_gly_dat <- fread("https://www.sciencebase.gov/catalog/file/get/5928a5dce4b016f7a93f8d7b?f=__disk__f7%2Ff0%2F8e%2Ff7f08ec5b316028dbe466c78e1159a196d1101db", skip = 1)

#Clean glyphosate and all other data and put together
msqa2013_gly <- msqa2013_gly_dat %>% 
  mutate(SAMPLE_DATE = as.Date(SAMPLE_DATE, format = "%m/%d/%Y"),
         PARM_NM = "Glyphosate, wf") %>% 
  rename("RESULT_CEN" = "GLYPHOSATE_ELISA_ng_L",
         "REMARK_CEN" = "GLYPHOSATE_RMK_ELISA") %>% 
  select(-c(SAMPLE_TIME_HHMM,MDL_ELISA_ng_L))

msqa2013_fin <- msqa2013_all_dat %>% 
  mutate(SAMPLE_DATE = as.Date(DATE2, format = "%m/%d/%Y")) %>% 
  select(SITE_NO, SHORT_NAME, SAMPLE_DATE, PARM_CD, PARM_NM, RESULT_CEN, REMARK_CEN) %>% 
  bind_rows(., msqa2013_gly) %>% 
  arrange(SAMPLE_DATE)

saveRDS(msqa2013_fin, "Agrochemical_Review/Sims/Data/midwest_pesticides_2013.rds")

# California SURF data from California Dept of Pesticide Regulation (CDPR) #########
ca_surf <- fread("ftp://transfer.cdpr.ca.gov/pub/outgoing/sw_data/SURF_Fusion18.csv")

ca_surf_clean <- ca_surf %>% 
  mutate(SAMP_DATE = ymd(SAMP_DATE))

saveRDS(ca_surf_clean, "Agrochemical_Review/Sims/Data/CA_SURF_clean.rds")
