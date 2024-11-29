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
species_code <- "KIPE" #four letter RAATD species code

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

#if n_users already equals 3, remove all variables to prevent further action
if(max(flag_IDs$n_users) >= 3){
  rm(list=ls())
  stop("Three users have already checked this species")
}

##------------------##
## User Input Start ##
##------------------##

# 3. Flag individuals with erroneous tracks

#define individual IDs flagged for removal
ind_flag <- c("king22_SG", "king24_SG")

#define corresponding device IDs flagged for removal - MUST BE SAME ORDER AS INDIVIDUAL IDs
dev_flag <- c("14354", "14360")

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

#if first user to flag for this species, create new dataframe and export
if(!exists("flag_IDs")){
  
  #select only combined IDs
  flagged <- flagged %>%
    select(combined_IDs) %>%
    mutate(combined_IDs = as.factor(combined_IDs))
  
  #add count and number of users columns 
  flagged <- flagged %>%
    mutate(count = 1,
           n_users = 1)
  
  #export
  saveRDS(flagged, file = "flags.RDS")
  
}

#if others have flagged already, append user flags
if(exists("flag_IDs")){
  
  #select only combined IDs
  flagged <- flagged %>% 
    select(combined_IDs) %>%
    mutate(combined_IDs = as.factor(combined_IDs))
  new_flags <- levels(flagged$combined_IDs)
  
  #add one to count if your flagged IDs are present already
  flag_IDs <- flag_IDs %>%
    mutate(count = if_else(combined_IDs %in% new_flags, count + 1, count))
  
  #extract IDs not previously flagged
  flag_IDs <- flag_IDs %>% 
    mutate(combined_IDs = as.factor(combined_IDs))
  old_flags <- levels(flag_IDs$combined_IDs)
  flagged <- flagged %>% 
    filter(!combined_IDs %in% old_flags) %>%
    mutate(count = 1,
           n_users = max(flag_IDs$n_users))
  flag_IDs <- flag_IDs %>% 
    bind_rows(flagged)
  
  #add one to number of users
  flag_IDs <- flag_IDs %>%
    mutate(n_users = n_users + 1)
  
  #export
  saveRDS(flag_IDs, file = "flags.RDS")
}


# 4. Remove flagged individuals from tracks and flag in metadata

#extract a vector of combined_IDs from flag_IDs where count >= 2
removal_IDs <- flag_IDs %>% 
  filter(count >= 2) %>%
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
  mutate(keepornot = if_else(combined_IDs %in% removal_IDs & is.na(keepornot), "discard", "keep"))

#remove combined_IDs column from tracks and metadata
tracks <- tracks %>% 
  select(-combined_IDs)
sp_meta <- sp_meta %>% 
  select(-combined_IDs)


# 5. Export dataframes

#export tracks 
saveRDS(tracks, file = paste0(species_code, "_ssm_qc.RDS"))

#replace keepornot column in meta with that from sp_meta
meta <- meta %>%
  filter(abbreviated_name != species_code) %>%
  mutate(keepornot = as.character()) %>%
  bind_rows(sp_meta)

#export updated metadata
saveRDS(meta, 
        file = "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
