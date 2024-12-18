#--------------------------------------------------------------
# Correct breeding stages based on visual checks and literature
#--------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
}


#----------------
# Adelie Penguins
#----------------

rm(list=ls())

this.species <- "ADPE"

#read in stages and tracks with stage trips
stages <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_stages.RDS"))
tracks <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_tracks_with_stage_trips.RDS"))

#append study codes from metadata
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% filter(abbreviated_name == this.species)


# 1. IAU_KG are all chick-rearing - see Machado-Gaye et al. (2024): doi.org/10.1007/s00227-024-04390-w

#correct iau_kgs to chick-rearing
iau_kg <- meta %>% filter(dataset_identifier == "IAU_KG")

#rename all stages to chick-rearing
iau_kg_stages <- stages %>% 
  filter(individual_id %in% iau_kg$individual_id) %>%
  mutate(stage = "chick-rearing")

#if multiple stages present, choose earliest start date and latest end date
iau_kg_stages <- iau_kg_stages %>%
  group_by(individual_id, device_id) %>%
  summarise(start = min(start), end = max(end), stage = "chick-rearing")


# 2. BAS_Signy_2008 are all correct already

#isolate bas_signy_2008 stages
bas_signy_2008 <- meta %>% filter(dataset_identifier == "BAS_Signy_2008")
bas_signy_2008_stages <- stages %>% 
  filter(individual_id %in% bas_signy_2008$individual_id) %>%
  filter(device_id %in% bas_signy_2008$device_id)


# 3. NIPR_Mamejima_2017 are all incubation - https://doi.org/10.1016/j.anbehav.2023.11.014 - already assigned as such

#isolate nipr_mamejima_2017 stages
nipr_mamejima_2017 <- meta %>% filter(dataset_identifier == "NIPR_Mamejima_2017")
nipr_mamejima_2017_stages <- stages %>% 
  filter(individual_id %in% nipr_mamejima_2017$individual_id) %>%
  filter(device_id %in% nipr_mamejima_2017$device_id)


# 4. USAMLR_WAP_2018 are all fledglings - https://doi.org/10.1098/rsbl.2020.0645 - temporal overlap with pre-moult

#isolate usamlr_wap_2018 stages
usamlr_wap_2018 <- meta %>% filter(dataset_identifier == "USAMLR_WAP_2018")
usamlr_wap_2018_stages <- stages %>% 
  filter(individual_id %in% usamlr_wap_2018$individual_id) %>%
  filter(device_id %in% usamlr_wap_2018$device_id)

#rename all stages to pre-moult
usamlr_wap_2018_stages <- usamlr_wap_2018_stages %>%
  mutate(stage = "pre-moult")


# 5. ADPE_CEBC_PG have breeding info in metadata - all are incubation or chick-rearing

#isolate adpe_cebc_pg stages
adpe_cebc_pg <- meta %>% filter(dataset_identifier == "ADPE_CEBC_PG")
adpe_cebc_pg_stages <- stages %>% 
  filter(individual_id %in% adpe_cebc_pg$individual_id) %>%
  filter(device_id %in% adpe_cebc_pg$device_id)

#isolate adpe_cebc_pg tracks
adpe_cebc_pg_tracks <- tracks %>% 
  filter(individual_id %in% adpe_cebc_pg$individual_id) %>%
  filter(device_id %in% adpe_cebc_pg$device_id)

#new chick-rearing IDs
chickrearing <- c("B12_AL169_CE_23−24_S1_2", "C1−CR AS88_S1_1", "C2−CR AS115_S1_1")

#new incubation IDs 
incubation <- c("A3_AL35_22−23_S1_1", "A5_AL40_22−23_S1_1")

#correct stages in tracks
adpe_cebc_pg_tracks <- adpe_cebc_pg_tracks %>%
  mutate(stage = case_when(
    trip %in% chickrearing ~ "chick-rearing",
    trip %in% incubation ~ "incubation",
    TRUE ~ stage
  ))

#reformat stages from tracks
adpe_cebc_pg_stages <- adpe_cebc_pg_tracks %>%
  group_by(individual_id, device_id, stage) %>%
  summarise(start = min(date), end = max(date))


# 6. NPOLAR_KG are all non-breeders (but some incubated first) - https://doi.org/10.1002/ecs2.4090

#isolate npolar_kg stages
npolar_kg <- meta %>% filter(dataset_identifier == "NPOLAR_KG")
npolar_kg_stages <- stages %>% 
  filter(individual_id %in% npolar_kg$individual_id) %>%
  filter(device_id %in% npolar_kg$device_id)

#isolate tracks for npolar_kg
npolar_kg_tracks <- tracks %>% 
  filter(individual_id %in% npolar_kg$individual_id) %>%
  filter(device_id %in% npolar_kg$device_id)

#new chick-rearing trip IDs
chickrearing <- c("153797_5", "153797_6", "153797_7", "153816_1", "153818_0",
                  "153818_1", "153819_3", "153819_4", "153819_5", "153823_1",
                  "153823_2", "153823_5", "153823_6", "153823_7", "153823_9",
                  "153801_5", "153803_2", "153803_4", "153803_5", "153811_7",
                  "153811_8", "153811_9", "153820_6", "153824_7", "153824_8",
                  "153824_9", "153825_4")

#new incubation trip IDs
incubation <- c("153807_1", "153808_1", "153813_2", "153813_4", "153814_1")

#new pre-moult trip IDs
premoult <- c("153797_11", "153798_3", "153802_4", "153807_6", "153808_14", "153809_6",
              "153810_3", "153813_6", "153814_4", "153816_2", "153818_11", "153819_6",
              "153822_9", "153823_18", "153826_2", "153799_7", "153800_7", "153801_11",
              "153803_17", "153804_6", "153805_13", "153806_4", "153812_18", "153815_4",
              "153824_14", "153825_10")

#correct stages in tracks
npolar_kg_tracks <- npolar_kg_tracks %>%
  mutate(stage = case_when(
    trip %in% chickrearing ~ "chick-rearing",
    trip %in% incubation ~ "incubation",
    trip %in% premoult ~ "pre-moult",
    TRUE ~ stage
  ))

#create stage dates from tracks
npolar_kg_stages <- npolar_kg_tracks %>%
  group_by(individual_id, device_id, stage) %>%
  summarise(start = min(date), end = max(date))


# 7. Join all together

#rejoin all stages
adpe_stages <- bind_rows(adpe_cebc_pg_stages, npolar_kg_stages, usamlr_wap_2018_stages,
                         bas_signy_2008_stages, iau_kg_stages, nipr_mamejima_2017_stages)

#export
saveRDS(adpe_stages, paste0("code/stage_splitting/stage_dates/", this.species, "_stages_corrected.RDS"))

#----------------
# Macaroni Penguins
#----------------

rm(list=ls())

this.species <- "MAPE"

#read in stages and tracks with stage trips
stages <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_stages.RDS"))
tracks <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_tracks_with_stage_trips.RDS"))

#append study codes from metadata
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% filter(abbreviated_name == this.species)


# 1. NMU_Marion_2018 breeding stages from seabirdtracking dataframe

#isolate nmu_marion_2018 stages
nmu_marion_2018 <- meta %>% filter(dataset_identifier == "NMU_Marion_2018")
nmu_marion_2018_stages <- stages %>% 
  filter(individual_id %in% nmu_marion_2018$individual_id) %>%
  filter(device_id %in% nmu_marion_2018$device_id)

#read in seabird tracking data
seabird <- read.csv(list.files(path = "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/MAPE/NMU_Marion_2018/Data", 
                               pattern = ".csv", full.names = T))

#format for later use
seabird <- seabird %>%
  rename(individual_id = bird_id,
         device_id = track_id) %>%
  mutate(date = paste0(date_gmt, time_gmt)) %>%
  mutate(date = as_datetime(date, format = "%d/%m/%Y %H:%M:%S"))

#find stage dates according to seabird tracking
seabird <- seabird %>%
  group_by(individual_id, device_id, breed_stage) %>%
  summarise(start = min(date), end = max(date)) %>%
  rename(stage = breed_stage)

#rename stages to RAATD names
seabird <- seabird %>%
  mutate(stage = case_when(
    stage == "incubation" ~ "incubation",
    stage == "brood-guard" ~ "early chick-rearing"
  ))

#bank these stages
nmu_marion_2018_stages <- seabird %>%
  mutate(individual_id = as.character(individual_id))


# 2. NMU_Marion_2015 are all early chick-rearing (from seabird tracking data info)

#isolate nmu_marion_2015 stages
nmu_marion_2015 <- meta %>% filter(dataset_identifier == "NMU_Marion_2015")
nmu_marion_2015_stages <- stages %>% 
  filter(individual_id %in% nmu_marion_2015$individual_id) %>%
  filter(device_id %in% nmu_marion_2015$device_id)

#rename all stages to early chick-rearing
nmu_marion_2015_stages <- nmu_marion_2015_stages %>%
  mutate(stage = "early chick-rearing")


# 3. BAS_SG_2011 individual_ids are coded by stage

#isolate bas_sg_2011 stages
bas_sg_2011 <- meta %>% filter(dataset_identifier == "BAS_SG_2011")
bas_sg_2011_stages <- stages %>% 
  filter(individual_id %in% bas_sg_2011$individual_id) %>%
  filter(device_id %in% bas_sg_2011$device_id)

#isolate bas_sg_2011 tracks
bas_sg_2011_tracks <- tracks %>% 
  filter(individual_id %in% bas_sg_2011$individual_id) %>%
  filter(device_id %in% bas_sg_2011$device_id)

#new early chick-rearing trip IDs
early <- c("b1a_1", "b1a_2", "b1b_1", "b1c_1", "b1d_1", "b1e_1", "c2a_2", "c3a_2", "c4a_1", "c4a_2",
           "c4b_1", "c4b_2", "c4b_3", "c4b_4", "c4c_1", "c4c_2", "c4c_3", "c4c_4", "i2b_0",
           "Obs140312_204206_Tag29330_1", "Obs140312_204557_Tag29318_1")

#new incubation trip IDs
incubation <- c("i1a_1", "i1b_1", "i1c_0")

#reassign trip stages
bas_sg_2011_tracks <- bas_sg_2011_tracks %>%
  mutate(stage = case_when(
    trip %in% early ~ "early chick-rearing",
    trip %in% incubation ~ "incubation",
    TRUE ~ stage
  ))

#make new stage dates
bas_sg_2011_stages <- bas_sg_2011_tracks %>%
  group_by(individual_id, device_id, stage) %>%
  summarise(start = min(date), end = max(date))


# 4. BAS_SG_2021 are all late chick-rearing (already assigned as that anyway)

#isolate bas_sg_2021 stages
bas_sg_2021 <- meta %>% filter(dataset_identifier == "BAS_SG_2021")
bas_sg_2021_stages <- stages %>% 
  filter(individual_id %in% bas_sg_2021$individual_id) %>%
  filter(device_id %in% bas_sg_2021$device_id)


# 5. Join all back together

#join all stages back together
mape_stages <- bind_rows(nmu_marion_2018_stages, nmu_marion_2015_stages, bas_sg_2011_stages, bas_sg_2021_stages)

#export corrected stages
saveRDS(mape_stages, paste0("code/stage_splitting/stage_dates/", this.species, "_stages_corrected.RDS"))


#----------------
# Emperor Penguins
#----------------

rm(list=ls())

this.species <- "EMPE"

#read in stages and tracks with stage trips
stages <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_stages.RDS"))
tracks <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_tracks_with_stage_trips.RDS"))

#append study codes from metadata
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% filter(abbreviated_name == this.species)


# 1. AWI_Atka_2019 are all fledglings - https://doi.org/10.1098/rsos.211708 (so post-breeding)

#isolate awi_atka_2019 stages
awi_atka_2019 <- meta %>% filter(dataset_identifier == "AWI_Atka_2019")
awi_atka_2019_stages <- stages %>% 
  filter(individual_id %in% awi_atka_2019$individual_id) %>%
  filter(device_id %in% awi_atka_2019$device_id)

#rename all stages to post-breeding
awi_atka_2019_stages <- awi_atka_2019_stages %>%
  mutate(stage = "post-breeding")


# 2. CEBC_PG_2013 are all fledglings - https://doi.org/10.3354/meps12831 (so post-breeding)

#isolate cebc_pg_2013 stages
cebc_pg_2013 <- meta %>% filter(dataset_identifier == "CEBC_PG_2013")
cebc_pg_2013_stages <- stages %>% 
  filter(individual_id %in% cebc_pg_2013$individual_id) %>%
  filter(device_id %in% cebc_pg_2013$device_id)

#rename all stages to post-breeding
cebc_pg_2013_stages <- cebc_pg_2013_stages %>%
  mutate(stage = "post-breeding")


# 3. CEBC_PG_1996 are all breeders - require visual state reassignment

#isolate cebc_pg_1996 stages
cebc_pg_1996 <- meta %>% filter(dataset_identifier == "CEBC_PG_1996")
cebc_pg_1996_stages <- stages %>% 
  filter(individual_id %in% cebc_pg_1996$individual_id) %>%
  filter(device_id %in% cebc_pg_1996$device_id)

#trips assigned as chick-rearing that don't return to the colony for months
nonbreeders <- c("26572c", "26574ab")

#reassign states for nonbreeders to post-breeding
cebc_pg_1996_stages <- cebc_pg_1996_stages %>%
  mutate(stage = case_when(
    individual_id %in% nonbreeders & stage == "chick-rearing" ~ "post-breeding",
    TRUE ~ stage
  ))


# 4. BAS_Rothschild_2015 breeding status unknown - require visual state reassignment

#isolate bas_rothschild_2015 stages
bas_rothschild_2015 <- meta %>% filter(dataset_identifier == "BAS_Rothschild_2015")
bas_rothschild_2015_stages <- stages %>% 
  filter(individual_id %in% bas_rothschild_2015$individual_id) %>%
  filter(device_id %in% bas_rothschild_2015$device_id)

#currently all assigned as chick-rearing but many trips don't return to colony
#list trips where individual does not return to the colony
postbreeders <- c("123005_3", "123006_2", "130997_1", "130999_5", "131000_2", "131001_1",
                  "131004_1", "131005_1", "131006_1", "131429_3", "131430_1", "131431_2",
                  "131432_1", "131433_8", "150624_1", "150626_3", "150630_2", "150631_2",
                  "98530_1", "98531_1", "98532_1", "98533_3", "98534_5", "98535_4",
                  "98537_2", "98539_1")

#isolate tracks for bas_rothschild_2015
bas_rothschild_2015_tracks <- tracks %>% 
  filter(individual_id %in% bas_rothschild_2015$individual_id) %>%
  filter(device_id %in% bas_rothschild_2015$device_id)

#reassign states for postbreeders to post-breeding
bas_rothschild_2015_tracks <- bas_rothschild_2015_tracks %>%
  mutate(stage = case_when(
    trip %in% postbreeders ~ "post-breeding",
    TRUE ~ stage
  ))

#convert updated trip info into breeding stages
bas_rothschild_2015_stages <- bas_rothschild_2015_tracks %>% 
  group_by(individual_id, device_id, stage) %>%
  summarise(start = first(date), end = last(date))


# 5. Join all back together

#join all stages back together
empe_stages <- bind_rows(awi_atka_2019_stages, cebc_pg_2013_stages, cebc_pg_1996_stages, bas_rothschild_2015_stages)

#export corrected stages
saveRDS(empe_stages, paste0("code/stage_splitting/stage_dates/", this.species, "_stages_corrected.RDS"))


#----------------
# King Penguins
#----------------

rm(list=ls())

this.species <- "KIPE"

#read in stages and tracks with stage trips
stages <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_stages.RDS"))
tracks <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_tracks_with_stage_trips.RDS"))

#append study codes from metadata
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% filter(abbreviated_name == this.species)


# 1. SAERI_FK are all chick-rearing - https://doi.org/10.1007/s00227-014-2561-0 - already assigned as such

#isolate saeri_fk stages
saeri_fk <- meta %>% filter(dataset_identifier == "SAERI_FK")
saeri_fk_stages <- stages %>% 
  filter(individual_id %in% saeri_fk$individual_id) %>%
  filter(device_id %in% saeri_fk$device_id)


# 2. SG_StAnd_2017 are all chick-rearing - https://doi.org/10.3389/fmars.2021.745200

#isolate sg_stand_2017 stages
sg_stand_2017 <- meta %>% filter(dataset_identifier == "SG_StAnd_2017")
sg_stand_2017_stages <- stages %>% 
  filter(individual_id %in% sg_stand_2017$individual_id) %>%
  filter(device_id %in% sg_stand_2017$device_id)

#change all stages to chick-rearing
sg_stand_2017_stages <- sg_stand_2017_stages %>%
  mutate(stage = "chick-rearing")


# 3. CNRS_Crozet_2012_2016 have no info available - visual reassignment

#isolate cnrs_crozet_2012_2016 stages
cnrs_crozet_2012_2016 <- meta %>% filter(dataset_identifier == "CNRS_Crozet_2012_2016")
cnrs_crozet_2012_2016_stages <- stages %>% 
  filter(individual_id %in% cnrs_crozet_2012_2016$individual_id) %>%
  filter(device_id %in% cnrs_crozet_2012_2016$device_id)

#isolate tracks for cnrs_crozet_2012_2016
tracks <- tracks %>% mutate(trip = as.character(trip))
cnrs_crozet_2012_2016_tracks <- tracks %>% 
  filter(individual_id %in% cnrs_crozet_2012_2016$individual_id) %>%
  filter(device_id %in% cnrs_crozet_2012_2016$device_id)

#a few chick-rearing trips are identified as incubation
chickrearing <- unique(cnrs_crozet_2012_2016_tracks$trip)[c(3, 4, 10, 26, 30, 35, 37, 38, 36)]

#a couple of incubation trips are identified as chick-rearing
incubation <- unique(cnrs_crozet_2012_2016_tracks$trip)[c(28, 34, 40)]

#reassign states for chickrearing to chick-rearing
cnrs_crozet_2012_2016_tracks <- cnrs_crozet_2012_2016_tracks %>%
  mutate(stage = ifelse(trip %in% chickrearing, "chick-rearing", stage)) %>%
  mutate(stage = ifelse(trip %in% incubation, "incubation", stage))

#convert updated trip info into breeding stages
cnrs_crozet_2012_2016_stages <- cnrs_crozet_2012_2016_tracks %>% 
  group_by(individual_id, device_id, stage) %>%
  summarise(start = min(date), end = max(date))


# 4. ART_Tierra_GPS and ART_Tierra_PTT are mostly correct but need some corrections - https://doi.org/10.1016/j.gecco.2021.e01669

#isolate art_tierra_gps and art_tierra_ptt stages
art_tierra <- meta %>% filter(dataset_identifier %in% c("ART_Tierra_GPS", "ART_Tierra_PTT"))
art_tierra_stages <- stages %>% 
  filter(individual_id %in% art_tierra$individual_id) %>%
  filter(device_id %in% art_tierra$device_id)

#pre-breeding tracks need changing to incubation
art_tierra_stages <- art_tierra_stages %>%
  mutate(stage = ifelse(stage == "Pre-breeding", "incubation", stage)) %>%
  group_by(individual_id, device_id, stage) %>%
  summarise(start = min(start), end = max(end)) %>%
  ungroup()

#two post-breeding tracks need correction
postbreeders <- c("146579b_15", "146579c_15")
art_tierra_stages <- art_tierra_stages %>%
  mutate(stage = ifelse(individual_id %in% postbreeders, "post-breeding", stage))


# 5. Join all back together

#join all stages back together
kipe_stages <- bind_rows(saeri_fk_stages, sg_stand_2017_stages, cnrs_crozet_2012_2016_stages, art_tierra_stages)

#export corrected stages
saveRDS(kipe_stages, paste0("code/stage_splitting/stage_dates/", this.species, "_stages_corrected.RDS"))


#----------------
# Chinstrap Penguins
#----------------

rm(list=ls())

this.species <- "CHPE"

#read in stages and tracks with stage trips
stages <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_stages.RDS"))
tracks <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_tracks_with_stage_trips.RDS"))

#append study codes from metadata
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% filter(abbreviated_name == this.species)


# 1. USAMLR_ADB_CHPE have breeding stages in seabirdtracking data

#isolate usamlr_adb_chpe stages
usamlr_adb_chpe <- meta %>% filter(dataset_identifier == "USAMLR_ADB_CHPE")
usamlr_adb_chpe_stages <- stages %>% 
  filter(individual_id %in% usamlr_adb_chpe$individual_id) %>%
  filter(device_id %in% usamlr_adb_chpe$device_id)

#read in seabird tracking data
seabird <- read.csv(list.files(path = "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/CHPE/USAMLR_ADB_CHPE/Data", 
                                       pattern = ".csv", full.names = T))

#format for later use
seabird <- seabird %>%
  rename(individual_id = bird_id,
         device_id = track_id) %>%
  mutate(date = paste(date_gmt, time_gmt, sep = " ")) %>%
  mutate(date = as_datetime(date, format = "%Y-%m-%d %H:%M:%S"))

#find stage dates according to seabird tracking
seabird <- seabird %>%
  group_by(individual_id, device_id, breed_stage) %>%
  summarise(start = min(date), end = max(date)) %>%
  rename(stage = breed_stage)

#rename stages to RAATD names
seabird <- seabird %>%
  mutate(stage = case_when(
    stage == "creche" ~ "chick-rearing",
    stage == "brood-guard" ~ "chick-rearing",
    stage == "winter" ~ "post-moult"
  ))

#filter to only remaining IDs
seabird <- seabird %>% filter(individual_id %in% usamlr_adb_chpe_stages$individual_id)

#bank these stages
usamlr_adb_chpe_stages <- seabird


# 2. USAMLR_CAS_CHPE have breeding stages in seabirdtracking data

#isolate usamlr_cas_chpe stages
usamlr_cas_chpe <- meta %>% filter(dataset_identifier == "USAMLR_CAS_CHPE")
usamlr_cas_chpe_stages <- stages %>% 
  filter(individual_id %in% usamlr_cas_chpe$individual_id) %>%
  filter(device_id %in% usamlr_cas_chpe$device_id)

#read in seabird tracking data
seabird <- read.csv(list.files(path = "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/CHPE/USAMLR_CAS_CHPE/Data", 
                                       pattern = ".csv", full.names = T))

#format for later use
seabird <- seabird %>%
  rename(individual_id = bird_id,
         device_id = track_id) %>%
  mutate(date = paste(date_gmt, time_gmt, sep = " ")) %>%
  mutate(date = as_datetime(date, format = "%d/%m/%Y %H:%M:%S"))

#find stage dates according to seabird tracking
seabird <- seabird %>%
  group_by(individual_id, device_id, breed_stage) %>%
  summarise(start = min(date), end = max(date)) %>%
  rename(stage = breed_stage)

#rename stages to RAATD names
seabird <- seabird %>%
  mutate(stage = case_when(
    stage == "creche" ~ "chick-rearing",
    stage == "brood-guard" ~ "chick-rearing",
    stage == "winter" ~ "post-moult",
    stage == "incubation" ~ "incubation"
  ))

#filter to only remaining IDs
seabird <- seabird %>%
  filter(individual_id %in% usamlr_cas_chpe_stages$individual_id)

#bank these stages
usamlr_cas_chpe_stages <- seabird


# 3. USAMLR_WAP_2017 need manual reassignment

#isolate usamlr_wap_2017 stages
usamlr_wap_2017 <- meta %>% filter(dataset_identifier == "USAMLR_WAP_2017")
usamlr_wap_2017_stages <- stages %>% 
  filter(individual_id %in% usamlr_wap_2017$individual_id) %>%
  filter(device_id %in% usamlr_wap_2017$device_id)

#isolate tracks for these individuals
usamlr_wap_2017_tracks <- tracks %>% filter(individual_id %in% usamlr_wap_2017$individual_id)

#all post-moult trip IDs
postmoult <- c("165206_7", "165211_15", "165212_1", "165213_7", "165271_1", "165273_1",
               "165151_4", "165153_7", "165154_2", "165156_14", "165158_9", "165159_1",
               "165160_2", "165161_2", "165162_1", "165164_11", "165257_1", "165188_12", 
               "165185_16", "165182_14", "165181_3", "165155_1", "165207_8", "165164_11",
               "165156_14", "165152_8")

#new incubation trip 
incubation <- c("165156_1")

#IDs that are mostly chick-rearing
chickrearing <- c("165182", "165185", "165188", "165152", "165153", "165154",
                  "165156", "165158", "165159", "165160", "165161", "165164",
                  "165206", "165207", "165211", "165213")

#reassign states for post-moult tracks
usamlr_wap_2017_tracks <- usamlr_wap_2017_tracks %>%
  mutate(stage = case_when(
    individual_id %in% chickrearing ~ "chick-rearing",
    TRUE ~ stage
  )) %>%
  mutate(stage = case_when(
    trip %in% postmoult ~ "post-moult",
    trip %in% incubation ~ "incubation",
    TRUE ~ stage
  ))

#convert updated trip info into breeding stages
usamlr_wap_2017_stages <- usamlr_wap_2017_tracks %>%
  group_by(individual_id, device_id, stage) %>%
  summarise(start = min(date), end = max(date))


# 4. BAS_South_Orkney_2011_2016 needs manual reassignment - use Polar Data Center metadata to cross-validate

#isolate bas_south_orkney_2011_2016 stages
bas_south_orkney_2011_2016 <- meta %>% filter(dataset_identifier == "BAS_South_Orkney_2011_2016")
bas_south_orkney_2011_2016_stages <- stages %>% 
  filter(individual_id %in% bas_south_orkney_2011_2016$individual_id) %>%
  filter(device_id %in% bas_south_orkney_2011_2016$device_id)

#isolate tracks for these individuals
bas_south_orkney_2011_2016_tracks <- tracks %>% 
  filter(individual_id %in% bas_south_orkney_2011_2016$individual_id)

#new chick-rearing IDs - CONTINUE
chickrearing <- c("LIR210_2", "LIR302_2", "LIR305_2", "LIR307_2", "LIR310_0", "LIR310_1",
                  "LIR401_3", "LIR402_1", "LIR403_3", "LIR403_4", "LIR403_5", "LIR404_3",
                  "LIR405_0", "LIR406_1", "LIR406_2", "LIR406_3", "LIR406_4", "LIR407_1",
                  "LIR408_1", "LIR409_1", "LIR409_2", "LIR409_3", "LIR410_1", "LIR410_2",
                  "LIR410_3", "LIR501_1", "LIR501_2", "LIR502_1", "LIR502_2", "LIR504_1",
                  "LIR506_1", "LIR507_1", "LIR508_1", "MIR301_1", "MIR301_2", "MIR302_1",
                  "PIR202_1", "PIR202_2", "PIR203_1", "PIR203_2", "PIR204_1", "PIR204_2",
                  "PIR204_3", "PIR205_1", "PIR207_1", "PIR207_2", "PIR208_0", "PIR208_1",
                  "PIR208_2", "PIR209_2", "PIR210_1", "PIR210_2", "PIR211_1", "PIR211_2",
                  "PIR212_1", "PIR212_2", "SIR103_2", "SIR112_1", "SIR112_2", "SIR112_3")

#new incubation IDs - CONTINUE
incubation <- c("SIR108_1", "SIR108_2", "SIR117_2")

#reassign trips to correct stages
bas_south_orkney_2011_2016_tracks <- bas_south_orkney_2011_2016_tracks %>%
  mutate(stage = case_when(
    trip %in% chickrearing ~ "chick-rearing",
    trip %in% incubation ~ "incubation",
    TRUE ~ stage
  ))

#convert updated trip info into breeding stages
bas_south_orkney_2011_2016_stages <- bas_south_orkney_2011_2016_tracks %>%
  group_by(individual_id, device_id, stage) %>%
  summarise(start = min(date), end = max(date))


# 5. BAS_Signy_2011_2016 has breeding stages in original metadata

#isolate bas_signy_2011_2016 stages
bas_signy_2011_2016 <- meta %>% filter(dataset_identifier == "BAS_Signy_2011_2016")
bas_signy_2011_2016_stages <- stages %>% 
  filter(individual_id %in% bas_signy_2011_2016$individual_id) %>%
  filter(device_id %in% bas_signy_2011_2016$device_id)

#new chick-rearing IDs
chickrearing <- c("29322_11_01_14", "29615_11_01_14")

#new incubation IDs
incubation <- c("29603_20_12_13", "29606_11_01_14", "29624_24_12_13")

#reassign individuals to correct stages
bas_signy_2011_2016_stages <- bas_signy_2011_2016_stages %>%
  mutate(stage = case_when(
    individual_id %in% chickrearing ~ "chick-rearing",
    individual_id %in% incubation ~ "incubation",
    TRUE ~ stage
  ))


# 6. BAS_Signy_2008 is already correct 

#isolate bas_signy_2008 stages
bas_signy_2008 <- meta %>% filter(dataset_identifier == "BAS_Signy_2008")
bas_signy_2008_stages <- stages %>% 
  filter(individual_id %in% bas_signy_2008$individual_id) %>%
  filter(device_id %in% bas_signy_2008$device_id)


# 7. BAS_Signy_2000 should all be chick-rearing (from seabirdtracking data)

#isolate bas_signy_2000 stages
bas_signy_2000 <- meta %>% filter(dataset_identifier == "BAS_Signy_2000")
bas_signy_2000_stages <- stages %>% 
  filter(individual_id %in% bas_signy_2000$individual_id) %>%
  filter(device_id %in% bas_signy_2000$device_id)

#rename all stages as chick-rearing
bas_signy_2000_stages <- bas_signy_2000_stages %>%
  mutate(stage = "chick-rearing")

#for any duplicates, take the earliest start and the latest end date
bas_signy_2000_stages <- bas_signy_2000_stages %>%
  group_by(individual_id, device_id, stage) %>%
  summarise(start = min(start), end = max(end))


# 8. Join all stages together

#join all stages together
chpe_stages <- bind_rows(bas_south_orkney_2011_2016_stages, bas_signy_2011_2016_stages, 
                        bas_signy_2008_stages, bas_signy_2000_stages,
                        usamlr_wap_2017_stages, usamlr_adb_chpe_stages,
                        usamlr_cas_chpe_stages)

#export corrected stages
saveRDS(chpe_stages, paste0("code/stage_splitting/stage_dates/", this.species, "_stages_corrected.RDS"))

