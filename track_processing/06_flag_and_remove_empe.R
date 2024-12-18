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
species_code <- "EMPE" #four letter RAATD species code

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
ind_flag <- c("150629", "26574c", "26575o", "26576i", "4631p", "4633n")

#define corresponding device IDs flagged for removal - MUST BE SAME ORDER AS INDIVIDUAL IDs
dev_flag <- c("150629", "26574c", "26575o", "26576i", "4631p", "4633n")

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

#identify remaining IDs
remaining_IDs <- tracks %>%
  mutate(combined_IDs = as.character(combined_IDs)) %>%
  select(combined_IDs) %>%
  distinct() %>%
  pull(combined_IDs)

#change final keepornot column
sp_meta <- sp_meta %>% 
  mutate(keepornot = if_else(combined_IDs %in% remaining_IDs, "keep", "discard"))

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
