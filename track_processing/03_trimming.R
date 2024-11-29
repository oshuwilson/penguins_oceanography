#trimming 
rm(list=ls())

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(ggmap)
  library(CCAMLRGIS)
}

##-------------------------------------##
## User Input Only Required for Step 1 ##
##-------------------------------------##

# 1. Read in standardised tracks and relevant metadata

#set required variables
species_code <- "ADPE"
study_code <- "ADPE_CEBC_PG"

##----------------##
## User Input End ##
##----------------##

#set working directory to draw raw files from
setwd(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code))

#read in tracks
tracks <- readRDS("standardised_tracks.RDS")

#read in and filter metadata
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% filter(dataset_identifier == study_code & abbreviated_name == species_code)


# 2. Trim tracks using interactive function

# function will load plots in X11 windows 
# plots display time since deployment against distance from deployment site
# click the start and end point of the section of the track to KEEP
# 
# exclude sections where:
# - tags are turned on early
# - animals remain at deployment site for extended period at the start/end
# - location quality or frequency deteriorates towards end of PTT tracks
# - animals remain in a small area for extended period at end of track
#
# returns filtered dataset without trimmed locations

#source function
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/functions/track_trim.R")

#if time gaps of over 50 days exist, split into trips by assigning new device id
tracks <- tracks %>%
  group_by(individual_id, device_id) %>%
  arrange(individual_id, device_id, date) %>%
  mutate(lagtime = lag(date)) %>%
  mutate(timediff = difftime(date, lagtime, units = "days")) %>%
  mutate(trip = ifelse(timediff > 50, 1, 0)) %>%
  mutate(trip = ifelse(is.na(trip), 0, trip)) %>%
  mutate(trip = cumsum(trip)) %>%
  mutate(device_id = if_else(trip != 0, paste0(device_id, "_", trip), device_id)) %>%
  select(-lagtime, -timediff, -trip) %>%
  ungroup()

#add new device_ids to metadata
new_devs <- levels(as.factor(tracks$device_id))
new_devs <- subset(new_devs, !new_devs %in% levels(as.factor(meta$device_id)))

if(length(new_devs) > 1){
  new_meta <- tracks %>%
    filter(device_id %in% new_devs) %>%
    group_by(individual_id, device_id) %>%
    summarise(deployment_date = first(date),
              deployment_decimal_latitude = first(lat),
              deployment_decimal_longitude = first(lon)) %>%
    mutate(deployment_year = year(deployment_date),
           deployment_month = month(deployment_date),
           deployment_day = day(deployment_date),
           deployment_time = as_hms(deployment_date)) %>%
    mutate(deployment_time = as.character(deployment_time)) %>%
    select(-deployment_date) %>% 
    left_join(select(meta, -device_id, -deployment_year, -deployment_month,
                     -deployment_day, -deployment_time, -deployment_decimal_latitude,
                     -deployment_decimal_longitude), 
              by = "individual_id")
  meta <- bind_rows(meta, new_meta)
}

#remove individuals with fewer than one day of data
long_ids <- tracks %>% 
  mutate(date_only = as_date(date)) %>%
  group_by(individual_id, device_id) %>%
  summarise(days = n_distinct(date_only)) %>%
  filter(days > 1) %>%
  ungroup() %>%
  mutate_if(is.factor, as.character) %>%
  mutate_if(is.character, as.factor)

tracks <- tracks %>% 
  filter(individual_id %in% long_ids$individual_id, device_id %in% long_ids$device_id) %>%
  mutate_if(is.factor, as.character) %>%
  mutate_if(is.character, as.factor)

#run function
tracks <- track_trim(tracks, meta)

#visualise
tracks_terra <- vect(tracks, geom = c("lon", "lat"), crs = "epsg:4326")
tracks_terra <- project(tracks_terra, "epsg:6932")

coast <- vect(load_Coastline())
crop_coast <- crop(coast, ext(tracks_terra))

ggplot() + geom_spatvector(data = tracks_terra, size = 0.5) +
  geom_spatvector(data = crop_coast) +
  theme_minimal()


# 3. Export

#format df for export
tracks <- tracks %>%
  select(-location_to_keep)

#export
saveRDS(tracks, paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code, "/trimmed_tracks.RDS"))

#update metadata for individuals removed
tracks <- tracks %>%
  mutate(individual_id = as.character(individual_id),
         device_id = as.character(device_id)) %>%
  mutate(individual_id = as.factor(individual_id),
         device_id = as.factor(device_id))

discards <- meta %>% 
  filter(!(individual_id %in% tracks$individual_id) & !(device_id %in% tracks$device_id)) %>%
  mutate(keepornot = "discard")

meta <- meta %>% anti_join(discards, by = c("individual_id", "device_id")) %>%
  bind_rows(discards)

all_meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
all_meta <- all_meta %>% 
  filter(dataset_identifier != study_code | abbreviated_name != species_code) %>% 
  bind_rows(meta)
saveRDS(all_meta, "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
saveRDS(meta, "Data/metadata.RDS")
