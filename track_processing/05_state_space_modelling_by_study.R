#code to state-space model all filtered tracks for a species
#cleanup
rm(list=ls())

#load libraries
{
  library(tidyverse)
  library(data.table)
  library(sf)
  library(terra)
  library(tidyterra)
  library(CCAMLRGIS)
  library(aniMotum)
  library(patchwork)
}

##------------------##
## User Input Start ##
##------------------##

# 1. Define Key Variables
species_code <- "KIPE" #four letter RAATD species code
study_code <- "CNRS_Kerguelen"

##----------------##
## User Input End ##
##----------------##

# 2. Read in prefiltered tracks

#set working directory
setwd(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code))

#read in tracks
tracks <- readRDS("filtered_tracks.RDS")

#read in metadata for this study
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% filter(dataset_identifier == study_code)

#obtain max velocity (m/s) from reference table
vmax <- read_csv("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/speed_filters.csv",
                 show_col_types = F) %>%
  filter(species == species_code) %>%
  select(max_speed) %>%
  as.numeric()

#filter tracks to this study
tracks <- tracks %>% 
  filter(individual_id %in% meta$individual_id) %>% 
  filter(device_id %in% meta$device_id)


# 3. Split tracks where gaps exist longer than two days

#calculate time difference
tracks <- tracks %>% 
  group_by(individual_id, device_id) %>%
  arrange(individual_id, device_id, date) %>% 
  mutate(lag_date = lag(date)) %>%
  mutate(timediff = difftime(date, lag_date, units = "days"))

#if time difference is greater than two days, begin new trip ID
tracks <- tracks %>% 
  mutate(timediff = ifelse(is.na(timediff), 0, as.numeric(timediff))) %>%
  mutate(trip_id = cumsum(timediff > 2)) %>% 
  ungroup()

# 3. Argos State-Space-Modelling with aniMotum

#if tracks are from different deployments on the same individual or vice versa, fit_ssm will fail
#thus a merged device/individual id is used for the purpose of fit_ssm when available
tracks <- tracks %>% 
  mutate(id = as.factor(paste(individual_id, device_id, trip_id, sep = "_"))) 

#remove tracks with fewer than 10 locations
tracks <- tracks %>% 
  group_by(id) %>% 
  filter(n() > 9) %>% 
  ungroup()

#fit tracks 
fit_trax <- fit_ssm(tracks,
                    time.step = 1, #time sampling interval in hours
                    vmax = vmax, #speed filter in m/s
                    model="crw")

#reroute tracks around land
fit_trax <- route_path(fit_trax)

# 4. Individual SSM Validation Plots

#source function
source("~/OneDrive - University of Southampton/Documents/Chapter 02/code/functions/mini_qc_plot.R")

#quality control plots
qcplot(fit = fit_trax, 
       spp = species_code, 
       device_type = "Mixed",  
       step_duration = 1,
       study_code = study_code)


# 5. Run quality control

#grab predicted locations from GPS fit
ssm_tracks <- grab(fit_trax, what = "predicted")

#join device ids with individual ids
individuals <- tracks %>% 
  select(id, individual_id, device_id, trip_id) %>% 
  group_by(id) %>% 
  summarise(individual_id = first(individual_id),
            device_id = first(device_id),
            trip_id = first(trip_id))
ssm_tracks <- ssm_tracks %>% left_join(individuals)

#select important columns and standardise names
ssm_tracks <- ssm_tracks %>% 
  ungroup() %>%
  select(individual_id, device_id, date, lon, lat, x.se, y.se) %>%
  rename(datetime = date, lon_se_km = x.se, lat_se_km = y.se)

#list individuals to remove
rm_individuals <- c()

#remove individuals from tracks
ssm_tracks <- ssm_tracks %>% 
  filter(!individual_id %in% rm_individuals)

#change metadata flags
meta <- meta %>%
  mutate(keepornot = ifelse(individual_id %in% as.character(ssm_tracks$individual_id) &
                              device_id %in% as.character(ssm_tracks$device_id), "keep", "remove"))

# 6. Export

#read in existing ssm_tracks
ssm_tracks_existing <- readRDS(paste0(species_code, "_ssm_qc.RDS"))

#rename datetime
ssm_tracks <- ssm_tracks %>% 
  rename(date = datetime)

#bind rows together
new_ssm_tracks <- bind_rows(ssm_tracks_existing, ssm_tracks)

#export
saveRDS(new_ssm_tracks,
        file = paste0(species_code, "_ssm_qc.RDS"))

#read in all metadata
allmeta <- readRDS("~/OneDrive - University of Southampton/Documents/Chapter 02/data/metadata.RDS")

#format meta
meta <- meta %>%
  select(-id)

#update keepornot
allmeta <- allmeta %>%
  filter(dataset_identifier != study_code) %>%
  bind_rows(meta)

#export
saveRDS(allmeta,
        file = "~/OneDrive - University of Southampton/Documents/Chapter 02/data/metadata.RDS")

#read in chapter 2 tracks
chap2 <- readRDS(paste0("~/OneDrive - University of Southampton/Documents/Chapter 02/data/tracks/", species_code, "/all_tracks.RDS"))

#join new tracks to chapter 2 tracks
ssm_tracks <- ssm_tracks %>%
  mutate(individual_id = as.character(individual_id), 
        device_id = as.character(device_id),
        longitude_se = NA,
        latitude_se = NA)

chap2 <- rbind(chap2, ssm_tracks)

#export
saveRDS(chap2,
        file = paste0("~/OneDrive - University of Southampton/Documents/Chapter 02/data/tracks/", species_code, "/all_tracks.RDS"))
