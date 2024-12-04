#code to state-space model all filtered tracks for a species
#split trips if time gap longer than X days?????

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

##----------------##
## User Input End ##
##----------------##

# 2. Read in prefiltered tracks

#set working directory
setwd(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code))

#read in tracks
tracks <- readRDS("filtered_tracks.RDS")

#obtain max velocity (m/s) from reference table
vmax <- read_csv("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/speed_filters.csv",
                 show_col_types = F) %>%
  filter(species == species_code) %>%
  select(max_speed) %>%
  as.numeric()


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
                    model="rw")

#reroute tracks around land
fit_trax <- route_path(fit_trax)

# 4. Individual SSM Validation Plots

#source function
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/functions/qcplot.R")

#quality control plots
qcplot(fit = fit_trax, 
       spp = species_code, 
       device_type = "Mixed",  
       step_duration = 1)


# 5. Export SSM tracks

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

#export
saveRDS(ssm_tracks,
        file = paste0(species_code, "_state_space_modelled.RDS"))
