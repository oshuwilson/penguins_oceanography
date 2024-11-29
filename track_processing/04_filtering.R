#filtering

#cleanup
rm(list=ls())

{
  library(tidyverse)
  library(data.table)
  library(terra)
  library(tidyterra)
  library(patchwork)
}

##-------------------------------------##
## User Input Required in Step 1 and 3 ##
##-------------------------------------##

# 1. Read in standardised tracks and relevant metadata

#set required variables
species_code <- "ADPE" #the 4-letter RAATD code
study_code <- "ADPE_CEBC_PG" #Dataset identifier
device_type <- "GPS" #either GPS, PTT, or GLS

##----------------##
## User Input End ##
##----------------##

#set working directory to draw raw files from
setwd(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code))

#read in tracks
tracks <- readRDS("trimmed_tracks.RDS")

#read in and filter metadata
all_meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- all_meta %>% filter(dataset_identifier == study_code & abbreviated_name == species_code)


# 2. Run prefilter function

#source function
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/functions/prefilter.R")

#run function
tracks <- prefilter(tracks, device_type = device_type)


# 3. Validation

#if lc unknown assign as "Z"
err <- tracks %>% filter(lc == "" | is.na(lc)) %>% mutate(lc="Z")
tracks <- tracks %>% 
  anti_join(err, by = c("individual_id", "device_id", "date", "lon", "lat")) %>%
  bind_rows(err)

#validate lc levels against list of acceptable levels
lc_levels <- c("G", "GL", "3", "2", "1", "0", "A", "B", "Z")
lc_error <- tracks %>% filter(!lc %in% lc_levels) 

#if prints "Inspect LC" there are lc levels that won't work in state-space modelling
if(nrow(lc_error) == 0) {
  print("Proceed")
} else {
  print("Inspect LC")
}

#cleanup
rm(lc_levels, lc_error)

#create plots to identify tracks with noisy estimates or temporally irregular locations
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/functions/noisechecker.R")
noisechecker(tracks)

##------------------##
## User Input Start ##
##------------------##

#flag combined_ids with erroneous tracks
# - noisy relative to track length
# - temporally irregular location estimates
flag <- c()

##----------------##
## User Input End ##
##----------------##

#remove individuals with noisy tracks
tracks <- tracks %>% 
  mutate(combined_id = as.factor(paste(individual_id, device_id, sep = "_"))) %>%
  filter(!combined_id %in% flag)

#update keepornot column in metadata
meta <- meta %>% 
  mutate(combined_id = as.factor(paste(individual_id, device_id, sep = "_"))) %>%
  mutate(keepornot = ifelse(combined_id %in% flag & dataset_identifier == study_code, "discard", keepornot)) %>%
  select(-combined_id)


# 4. Export

#import existing filtered tracks for this species and merge if applicable
try({
  spp_tracks <- readRDS(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/filtered_tracks.RDS"))
  tracks <- bind_rows(spp_tracks, tracks)
  tracks <- tracks %>% distinct(.keep_all = T) #check that no repeats are included
})

#export all tracks for this species
saveRDS(tracks, file = paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/filtered_tracks.RDS"))

#export updated metadata with keepornot column
all_meta <- all_meta %>% 
  filter(dataset_identifier != study_code | abbreviated_name != species_code) %>%
  bind_rows(meta) 
saveRDS(all_meta, file = "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
saveRDS(meta, "Data/metadata.RDS")
