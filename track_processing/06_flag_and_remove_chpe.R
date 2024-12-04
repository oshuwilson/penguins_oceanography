#code to flag individuals during quality control
#remove tracks and flag metadata where 2 or more users have flagged the same individual

#cleanup
rm(list=ls())

#load libraries
{
  library(tidyverse)
}

##------------------##
## User Input Start ##
##------------------##

# 1. Define Key Variables
species_code <- "CHPE" #four letter RAATD species code

##----------------##
## User Input End ##
##----------------##

# 2. Read in dataframes

#set working directory
setwd(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code))

#load metadata for this species
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
sp_meta <- meta %>% filter(abbreviated_name == species_code)

#read in state-space-modelled tracks
tracks <- readRDS(paste0(species_code, "_state_space_modelled.RDS"))

#read in flagging data if already exists
try(flag_IDs <- readRDS("flags.RDS"))

#if flag IDs exist, clean environment to prevent further action
if(exists("flag_IDs")){
  rm(list=ls())
  stop("This species has already been checked")
}

##------------------##
## User Input Start ##
##------------------##

# 3. Flag individuals with erroneous tracks

#define individual IDs flagged for removal
ind_flag <- c("165150", "165157", "165163", "165180", "165183", "165186", "165189", "165205",
              "165208", "165209", "165210", "165214", "165255", "165256", "165258", "165270",
              "18", "47", "C008_23_01_08", "CHPE1", "CHPE1", 
              "CHPE101", "CHPE104", "CHPE105", "CHPE113",
              "CHPE116", "CHPE120", "CHPE121", "CHPE123",
              "CHPE124", "CHPE126", "CHPE127", "CHPE128",
              "CHPE132", "CHPE133", "CHPE138", "CHPE14",
              "CHPE142", "CHPE148", "CHPE150", "CHPE153",
              "CHPE173", "CHPE174", "CHPE180", "CHPE186",
              "CHPE189", "CHPE190", "CHPE195", "CHPE199",
              "CHPE2", "CHPE202", "CHPE208", "CHPE220",
              "CHPE23", "CHPE26", "CHPE27", "CHPE3",
              "CHPE30", "CHPE34", "CHPE36", "CHPE38",
              "CHPE4", "CHPE40", "CHPE44", "CHPE46",
              "CHPE51", "CHPE53", "CHPE69", "CHPE74", 
              "CHPE87", "F6", "J8", "K2",
              "LIR102", "M1", "M3")

#define corresponding device IDs flagged for removal - MUST BE SAME ORDER AS INDIVIDUAL IDs
dev_flag <- c("165150", "165157", "165163", "165180", "165183", "165186", "165189", "165205",
              "165208", "165209", "165210", "165214", "165255", "165256", "165258", "165270",
              "765_18_11002", "765_47_11015", "GPL008", "911_CHPE1_16674", "913_CHPE1_16823",
              "913_CHPE101_16826", "913_CHPE104_16829", "913_CHPE105_16830", "913_CHPE113_16839",
              "913_CHPE116_16842", "913_CHPE120_16847", "913_CHPE121_16848", "913_CHPE123_16850",
              "913_CHPE124_16851", "913_CHPE126_16853", "913_CHPE127_16854", "913_CHPE128_16855",
              "913_CHPE132_16860", "913_CHPE133_16861", "913_CHPE138_16866", "911_CHPE14_16679",
              "913_CHPE142_16871", "913_CHPE148_16877", "913_CHPE150_16880", "913_CHPE153_16883",
              "913_CHPE173_16905", "913_CHPE174_16906", "913_CHPE180_16913", "913_CHPE186_16919",
              "913_CHPE189_16922", "913_CHPE190_16924", "913_CHPE195_16929", "913_CHPE199_16933",
              "911_CHPE2_16685", "913_CHPE202_16938", "913_CHPE208_16944", "913_CHPE220_16958",
              "913_CHPE23_16968", "913_CHPE26_16971", "913_CHPE27_16972", "913_CHPE3_16975",
              "913_CHPE30_16976", "913_CHPE34_16980", "911_CHPE36_16703", "911_CHPE38_16705",
              "913_CHPE4_16986", "911_CHPE40_16708", "911_CHPE44_16712", "911_CHPE46_16714",
              "911_CHPE51_16720", "913_CHPE53_17001", "913_CHPE69_17018", "913_CHPE74_17024",
              "913_CHPE87_17038", "765_F6_11027", "765_J8_11032", "765_K2_11033",
              "LIR102", "765_M1_11040", "765_M3_11042")

##----------------##
## User Input End ##
##----------------##

#combine individual and device flag columns
flagged <- data.frame(individual = ind_flag, device = dev_flag)

#create merged column to ensure duplicated individuals or devices aren't also removed
flagged <- flagged %>% 
  mutate(combined_IDs = paste(individual, device, sep = "_"))

#check that combined_IDs match those in the metadata and identify any that don't
sp_meta <- sp_meta %>% 
  mutate(combined_IDs = as.factor(paste(individual_id, device_id, sep = "_")))

if(!all(flagged$combined_IDs %in% sp_meta$combined_IDs)){
  error <- flagged %>% 
    filter(!combined_IDs %in% sp_meta$combined_IDs) %>%
    select(combined_IDs) 
  print(paste0("The following combined IDs were not found in the metadata: ", error))
  print("Please check for spelling mistakes in ind_flag or dev_flag")
  print("Also check that individual and device IDs have been written in the same order")
  stop("Read above messages")
}

#create dataframe of flagged IDs and export
flagged <- flagged %>%
  select(combined_IDs) %>%
  mutate(combined_IDs = as.factor(combined_IDs))

#add count and number of users columns 
flagged <- flagged %>%
  mutate(count = 1,
         n_users = 1)

#export
saveRDS(flagged, file = "flags.RDS")


# 4. Remove flagged individuals from tracks and flag in metadata

#extract a vector of combined_IDs from flag_IDs where count >= 2
removal_IDs <- flagged %>%
  select(combined_IDs) %>%
  mutate(combined_IDs = as.character(combined_IDs)) %>%
  as.vector() %>%
  unlist()

#create combined_IDs column in tracks
tracks <- tracks %>% 
  mutate(combined_IDs = as.factor(paste(individual_id, device_id, sep = "_")))

#remove tracks from individuals that have been flagged
tracks <- tracks %>% 
  filter(!combined_IDs %in% removal_IDs)

#flag metadata for individuals that have been flagged
sp_meta <- sp_meta %>% 
  mutate(keepornot = if_else(!combined_IDs %in% removal_IDs & is.na(keepornot), "keep", "discard"))

#remove combined_IDs column from tracks and metadata
tracks <- tracks %>% 
  select(-combined_IDs)
sp_meta <- sp_meta %>% 
  select(-combined_IDs)


# 5. Export dataframes

#change datetime to date (to match RAATD 1.0)
tracks <- tracks %>% 
  rename(date = datetime)

#export tracks 
saveRDS(tracks, file = paste0(species_code, "_ssm_qc.RDS"))

#replace keepornot column in meta with that from sp_meta
meta <- meta %>%
  filter(abbreviated_name != species_code) %>%
  mutate(keepornot = as.character(keepornot)) %>%
  bind_rows(sp_meta)

#export updated metadata
saveRDS(meta, 
        file = "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
