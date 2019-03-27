#Data load and clean

require(tidyverse)

#My interventions data generated form interpretation of spreadsheets shared by Giulio
interventions <- read_csv("~/RemaisWork/Schisto/Stanford/Human_Parasitology_data/Village_interventions.csv") %>% 
  mutate(Intervention = factor(Intervention, levels = c("control", "net", "prawn", "veg_removal")),
         Region = factor(Region),
         Prawn_cat = factor(Prawn_cat, levels = c("No", "low", "high")))

dat1 <- read_csv("~/RemaisWork/Schisto/Stanford/Human_Parasitology_data/First_round_human_subject_analysis_APLS_POC_3.csv", na = c("DM", "NA"))
#data munging
  #sort(unique(dat1$P2_omega_total))

# Get means from different samples, rename variables, and reorder    
  dat1 <- dat1 %>% 
    group_by(Child_ID) %>% 
    mutate(year = 2016, School = as.character(School),
          # Replace >200 with 200 then coerce to numeric
           P2_omega_total = as.numeric(ifelse(P2_omega_total == ">200", 200, P2_omega_total)),
           #Number of haematobium samples
           haem_samps = length(na.omit(c(P1_omega_total, P2_omega_total))),
          #Mean haematobium eggs/10mL
           s_haem_mean2016 = mean(c(P1_omega_total, P2_omega_total)),  
          #Mean haematobium eggs/10mL NAs removed
           s_haem_mean2016_narm = mean(c(P1_omega_total, P2_omega_total), na.rm = TRUE), 
          #Number of samples contributing to mean
           s_haem_mean2016_w = length(which(!is.na(c(P1_omega_total, P2_omega_total)))),
          #WHO infection intensity labels
           s_haem_mean_cat2016 = cut(s_haem_mean2016_narm, 
                                     breaks = c(0, 1, 50, Inf), 
                                     labels = c("No", "Lo", "Hi"), right = FALSE), 
          #Mean mansoni epg from 1st kato katz
           s_mans_mean1 = mean(c(P1_kk1, P1_kk2)),           
          ##Mean mansoni epg from 2nd kato katz
           s_mans_mean2 = mean(c(P2_kk1, P2_kk2)),    
          #Mean mansoni epg from 1st/2nd kato katz
           s_mans_mean2016 = mean(c(P1_kk1, P1_kk2, P2_kk1, P2_kk2)), 
          #Mean mansoni EPG NAs removed
           s_mans_mean2016_narm = mean(c(P1_kk1, P1_kk2, P2_kk1, P2_kk2), na.rm = TRUE),  
          #Number of non-NA measurements
           s_mans_mean2016_w = length(which(!is.na(c(P1_kk1, P1_kk2, P2_kk1, P2_kk2)))),  
          #WHO infection levels
           s_mans_mean_cat2016 = cut(s_mans_mean2016_narm, 
                                     breaks = c(0,1,100,400,Inf), 
                                     labels = c("No", "Lo", "Md", "Hi"), right = FALSE),
           pzq = 0, 
           extra_pzq = 0) %>% 
    select(Child_ID, year, School, pzq, extra_pzq,
           s_haem1 = P1_omega_total,
           s_haem2 = P2_omega_total,
           haem_samps,
           s_haem_mean = s_haem_mean2016,
           s_haem_mean_narm = s_haem_mean2016_narm,
           s_haem_mean_w = s_haem_mean2016_w,
           s_haem_mean_cat = s_haem_mean_cat2016,
           s_mans1_1 = P1_kk1,
           s_mans1_2 = P1_kk2,
           s_mans2_1 = P2_kk1,
           s_mans2_2 = P2_kk2,
           s_mans_mean1,
           s_mans_mean2,
           s_mans_mean12 = s_mans_mean2016,
           s_mans_mean_narm = s_mans_mean2016_narm,
           s_mans_mean_w = s_mans_mean2016_w,
           s_mans_mean_cat = s_mans_mean_cat2016) %>% ungroup()
  
#Data check
  #make sure each subject only has one record:
    bad_id_dat1 <- dat1 %>% 
      group_by(Child_ID) %>% 
      summarise(nobs = n()) %>% 
      filter(nobs > 1) %>% 
      pull(Child_ID)
    
    length(bad_id_dat1)

dat2 <- read_csv("~/RemaisWork/Schisto/Stanford/Human_Parasitology_data/Second_round_human_subject_analysis_APLS_POC_2017_2.csv", na = c("DM", "SE", "NA"))

# Get means from different samples, rename variables, and reorder    
  dat2 <- dat2 %>% 
    group_by(Child_ID) %>% 
    mutate(year = 2017, School = as.character(School),
          # Replace >200 with 200 then coerce to numeric
           P4_omega_total = as.numeric(ifelse(P4_omega_total == ">200", 200, P4_omega_total)),
           #Number of haematobium samples
           haem_samps = length(na.omit(c(P3_omega_total, P4_omega_total))),
           #Mean haematobium eggs/10mL
           s_haem_mean2017 = mean(c(P3_omega_total, P4_omega_total)),  
           #Mean haematobium eggs/10mL NAs removed
           s_haem_mean2017_narm = mean(c(P3_omega_total, P4_omega_total), na.rm = TRUE), 
           #Number of samples contributing to mean
           s_haem_mean2017_w = length(which(!is.na(c(P3_omega_total, P4_omega_total)))),
          #WHO infection intensity labels
           s_haem_mean_cat2017 = cut(s_haem_mean2017_narm, 
                                     breaks = c(0, 1, 50, Inf), 
                                     labels = c("No", "Lo", "Hi"), right = FALSE), 
          #Mean mansoni epg from 1st kato katz
           s_mans_mean1 = mean(c(P3_kk1, P3_kk2)), 
          #Mean mansoni epg from 2nd kato katz
           s_mans_mean2 = mean(c(P4_kk1, P4_kk2)),      
          #Mean mansoni epg from 1st/2nd kato katz
           s_mans_mean2017 = mean(c(P3_kk1, P3_kk2, P4_kk1, P4_kk2)), 
          #Mean mansoni epg NAs removed
           s_mans_mean2017_narm = mean(c(P3_kk1, P3_kk2, P4_kk1, P4_kk2), na.rm = TRUE),  
          #Number of non-NA measurements
           s_mans_mean2017_w = length(which(!is.na(c(P3_kk1, P3_kk2, P4_kk1, P4_kk2)))),  
          #Number of samples contributing to mean
           s_mans_mean_cat2017 = cut(s_mans_mean2017_narm, 
                                     breaks = c(0,1,100,400, Inf), 
                                     labels = c("No", "Lo", "Md", "Hi"), right = FALSE),
           pzq = 1,
           extra_pzq = as.numeric(!is.na(PZQ_date_before_analyses)))  %>% 
    select(Child_ID, year, School, pzq, extra_pzq,
           s_haem1 = P3_omega_total,
           s_haem2 = P4_omega_total,
           haem_samps,
           s_haem_mean = s_haem_mean2017,
           s_haem_mean_narm = s_haem_mean2017_narm,
           s_haem_mean_w = s_haem_mean2017_w,
           s_haem_mean_cat = s_haem_mean_cat2017,
           s_mans1_1 = P3_kk1,
           s_mans1_2 = P3_kk2,
           s_mans2_1 = P4_kk1,
           s_mans2_2 = P4_kk2,
           s_mans_mean1,
           s_mans_mean2,
           s_mans_mean12 = s_mans_mean2017,
           s_mans_mean_narm = s_mans_mean2017_narm,
           s_mans_mean_w = s_mans_mean2017_w,
           s_mans_mean_cat = s_mans_mean_cat2017) %>% ungroup()
  
#Data check
  #make sure each subject only has one record:
    bad_id_dat2 <- dat2 %>% 
      group_by(Child_ID) %>% 
      summarise(nobs = n()) %>% 
      filter(nobs > 1) %>% 
      pull(Child_ID)
    
    length(bad_id_dat2)

dat3 <- read_csv("~/RemaisWork/Schisto/Stanford/Human_Parasitology_data/Third_round_human_subject_analysis_APLS_POC_2018_2.csv", na = c("DM", "SE", "NA"))
# Replace **_ numbers with just number and coerce to numeric
  dat3$P5_omega_total <- as.numeric(sub("**", "", dat3$P5_omega_total, fixed = TRUE))
  dat3$P6_omega_total <- sub("**", "*", dat3$P6_omega_total, fixed = TRUE)
  dat3$P6_omega_total <- as.numeric(sub("*", "", dat3$P6_omega_total, fixed = TRUE))

# Get means from different samples, rename variables, and reorder    
  dat3 <- dat3 %>% 
    group_by(Child_ID) %>% 
    mutate(year = 2018,
           School = strsplit(Child_ID, "/")[[1]][1],
           #Number of haematobium samples
           haem_samps = length(na.omit(c(P5_omega_total, P6_omega_total))),
           #Mean haematobium eggs/10mL
           s_haem_mean2018 = mean(c(P5_omega_total, P6_omega_total)),  
           #Mean haematobium eggs/10mL, NAs removed
           s_haem_mean2018_narm = mean(c(P5_omega_total, P6_omega_total), na.rm = TRUE), 
           #Number of samples contributing to mean
           s_haem_mean2018_w = length(which(!is.na(c(P5_omega_total, P6_omega_total)))),
           #WHO infection intensity labels
           s_haem_mean_cat2018 = cut(s_haem_mean2018_narm, 
                                     breaks = c(0, 1, 50, Inf), 
                                     labels = c("No", "Lo", "Hi"), right = FALSE), 
           #Mean mansoni epg from 1st kato katz
           s_mans_mean1 = mean(c(P5_kk1_epg, P5_kk2_epg)),   
           ##Mean mansoni epg from 2nd kato katz
           s_mans_mean2 = mean(c(P6_kk1_epg, P6_kk2_epg)),        
           #Mean mansoni epg from 1st/2nd kato katz
           s_mans_mean2018 = mean(c(P5_kk1_epg, P5_kk2_epg, P6_kk1_epg, P6_kk2_epg)), 
           #Mean mansoni epg from 1st/2nd kato katz, NAs removed
           s_mans_mean2018_narm = mean(c(P5_kk1_epg, P5_kk2_epg, P6_kk1_epg, P6_kk2_epg), na.rm = TRUE), 
           #Number of samples contributing to mean
           s_mans_mean2018_w = length(which(!is.na(c(P5_kk1_epg, P5_kk2_epg, P6_kk1_epg, P6_kk2_epg)))),
           #WHO infection categories
           s_mans_mean_cat2018 = cut(s_mans_mean2018_narm, 
                                     breaks = c(0,1,100,400, Inf), 
                                     labels = c("No", "Lo", "Md", "Hi"), right = FALSE),
           pzq = 1, 
           #Individuals who received additional praziquantel influencing their 2018 egg burden
           extra_pzq = ifelse(Child_ID %in% c("DT/1/008", "DT/1/021","DT/1/032",
                                              "MD/3/015", "MD/3/023",
                                              "MO/2/013", "MO/2/014", "MO/3/030"), 1, 0))  %>% 
    select(Child_ID, year, School, pzq, extra_pzq,
           s_haem1 = P5_omega_total,
           s_haem2 = P6_omega_total,
           haem_samps,
           s_haem_mean = s_haem_mean2018,
           s_haem_mean_narm = s_haem_mean2018_narm,
           s_haem_mean_w = s_haem_mean2018_w,
           s_haem_mean_cat = s_haem_mean_cat2018,
           s_mans1_1 = P5_kk1_epg,
           s_mans1_2 = P5_kk2_epg,
           s_mans2_1 = P6_kk1_epg,
           s_mans2_2 = P6_kk2_epg,
           s_mans_mean1,
           s_mans_mean2,
           s_mans_mean12 = s_mans_mean2018,
           s_mans_mean_narm = s_mans_mean2018_narm,
           s_mans_mean_w = s_mans_mean2018_w,
           s_mans_mean_cat = s_mans_mean_cat2018) %>% ungroup()
  
#Data check
  #make sure each subject only has one record:
    bad_id_dat3 <- dat3 %>% 
      group_by(Child_ID) %>% 
      summarise(nobs = n()) %>% 
      filter(nobs > 1) %>% 
      pull(Child_ID)
    
      bad_id_dat3
      
  #This person ended up with two records somehow, going to take the safe route: remove both records
    dat3 <- dat3 %>% filter(Child_ID != enquo(bad_id_dat3)[[2]])

#Full dataset (in long format) shared by Stanford, constructed by Isabel
stanford <- read_csv("~/RemaisWork/Schisto/Stanford/Human_Parasitology_data/human_data_2018Mat1518.csv") %>% 
  mutate(ID = factor(ID),    #Convert some variables to factors
         School = factor(School),
         Village = factor(Village),
         year_fac = factor(year),
         cluster = factor(cluster),
         Intervention.type = factor(Intervention.type),
         pzq = ifelse(year %in% c(2017, 2018), 1, 0))

#same code as above to remove individual with multiple observation in the same year
    stanford <- stanford %>% filter(ID != enquo(bad_id_dat3)[[2]])

#Jason's full dataset 
jason_dat <- read_csv("~/RemaisWork/Schisto/Stanford/Human_Parasitology_data/Jason_human_data_2018.csv") %>% 
    mutate(ID = factor(ID),    #Convert some variables to factors
           School = factor(School),
           Village = factor(Village),
           year_fac = factor(year),
           sex = factor(sex),
           cluster = factor(cluster),
           Intervention.type = factor(Intervention.type),
           pzq = ifelse(year %in% c(2017, 2018), 1, 0))


#Full dataset in long format from data reads/manipulations above
full_long <- rbind(dat1, dat2, dat3) %>% 
  full_join(interventions, by = c("School", "year")) %>% 
  #Join with stanford dataset to get cluster and sex
  full_join(stanford %>% select(ID, year, sex, cluster),    
            by = c("Child_ID" = "ID", "year" = "year")) %>% 
  mutate(School = factor(School),
         Child_ID = factor(Child_ID),
         s_haem_bin_narm = ifelse(s_haem_mean_narm >0, 1, 0),
         s_mans_bin_narm = ifelse(s_mans_mean_narm >0, 1, 0),
         year_fac = factor(year),
         Intervention_2018 = factor(Intervention_2018,
                                    levels = c("control", "prawn", "veg_removal")),
         year_2016 = ifelse(year == 2016, 1, 0),
         year_2017 = ifelse(year == 2017, 1, 0),
         year_2018 = ifelse(year == 2018, 1, 0),
         pzq = factor(pzq))

#Full community-level dataset in long format
full_long_comm <- full_long %>% #filter(extra_pzq != 1) %>%
  group_by(School, year) %>% 
  summarise(haem_mean = mean(c(s_haem1, s_haem2), na.rm = TRUE),
            haem_samps = length(na.omit(c(s_haem1, s_haem2))),
            haem_pos = length(which(s_haem_mean_narm > 0)), 
            haem_prev = haem_pos / length(na.omit(s_haem_mean_narm)),
            mans_mean = mean(c(s_mans1_1, s_mans1_2,
                               s_mans2_1, s_mans2_2), na.rm = TRUE),
            mans_samps = length(na.omit(c(s_mans1_1, s_mans1_2,
                                          s_mans2_1, s_mans2_2))),
            mans_pos = length(which(s_mans_mean_narm > 0)),
            mans_prev = mans_pos / length(na.omit(s_mans_mean_narm)),
            extra_pzq = mean(extra_pzq)) %>% 
  inner_join(interventions, by = c("School", "year")) %>% 
  group_by(School) %>% 
  mutate(vil_haem_mean = mean(haem_mean),
         vil_mans_mean = mean(mans_mean)) %>% 
  ungroup() %>% 
  mutate(vil_haem_diff = haem_mean - vil_haem_mean,
         vil_mans_diff = mans_mean - vil_mans_mean) %>% 
  arrange(year, School)

#Dataset with extra_pzq kids and kids with fewer samples removed
full_long_comm2 <- full_long %>% filter(extra_pzq != 1) %>%
  group_by(School, year) %>% 
  summarise(haem_mean = mean(c(s_haem1, s_haem2), na.rm = TRUE),
            haem_samps = length(na.omit(c(s_haem1, s_haem2))),
            haem_pos = length(which(s_haem_mean_narm > 0)), 
            haem_prev = haem_pos / length(na.omit(s_haem_mean_narm)),
            mans_mean = mean(c(s_mans1_1, s_mans1_2,
                               s_mans2_1, s_mans2_2), na.rm = TRUE),
            mans_samps = length(na.omit(c(s_mans1_1, s_mans1_2,
                                          s_mans2_1, s_mans2_2))),
            mans_pos = length(which(s_mans_mean_narm > 0)),
            mans_prev = mans_pos / length(na.omit(s_mans_mean_narm)),
            extra_pzq = mean(extra_pzq)) %>% 
  inner_join(interventions, by = c("School", "year")) %>% 
  group_by(School) %>% 
  mutate(vil_haem_mean = mean(haem_mean),
         vil_mans_mean = mean(mans_mean)) %>% 
  ungroup() %>% 
  mutate(vil_haem_diff = haem_mean - vil_haem_mean,
         vil_mans_diff = mans_mean - vil_mans_mean) %>% 
  arrange(year, School)

#Full dataset in wide format
full_wide <- dat1 %>% 
  full_join(dat2 %>% filter(extra_pzq != 1), 
            by = c("Child_ID", "School"), 
            suffix = c("_2016", "_2017")) %>% 
  full_join(dat3, by = c("Child_ID", "School")) %>% 
  full_join(interventions %>% filter(year == 2018), by = c("School"))

