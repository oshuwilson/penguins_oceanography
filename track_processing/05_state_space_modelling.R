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
species_code <- "CHPE" #four letter RAATD species code

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


# 3. Argos State-Space-Modelling with aniMotum

#if tracks are from different deployments on the same individual or vice versa, fit_ssm will fail
#thus a merged device/individual id is used for the purpose of fit_ssm when available
tracks <- tracks %>% 
  mutate(id = as.factor(paste0(individual_id, "_", device_id))) 

#split tracks by device_type before state-space modelling (different time steps for each)
GPS <- tracks %>% 
  filter(device_type == "GPS")

PTT <- tracks %>% 
  filter(device_type == "PTT")

GLS <- tracks %>%
  filter(device_type == "GLS")

#fit GPS tracks and reroute around land
if(nrow(GPS) > 0){
  fit_GPS <- fit_ssm(GPS,
                     time.step = 1, #time sampling interval in hours
                     vmax = vmax, #speed filter in m/s
                     model="crw")
  fit_GPS <- route_path(fit_GPS)
}

#repeat for PTT
if(nrow(PTT) > 0){
  fit_PTT <- fit_ssm(PTT,
                     time.step = 2,
                     vmax = vmax,
                     model = "crw")
  fit_PTT <- route_path(fit_PTT)
}

#repeat for GLS
if(nrow(GLS) > 0){
  fit_GLS <- fit_ssm(GLS,
                     time.step = 12,
                     vmax = vmax,
                     model = "crw")
  fit_GLS <- route_path(fit_GLS)
}


# 4. Individual SSM Validation Plots

#source function
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/functions/qcplot.R")

#GPS plots
qcplot(fit = fit_GPS, 
       spp = species_code, 
       device_type = "GPS", 
       step_duration = 1)

#PTT plots
qcplot(fit = fit_PTT, 
       spp = species_code, 
       device_type = "PTT", 
       step_duration = 2)

#GLS plots
qcplot(fit = fit_GLS,
       spp = species_code,
       device_type = "GLS",
       step_duration = 12)


# 5. Export SSM tracks

#grab predicted locations from GPS fit
GPS_ssm_tracks <- grab(fit_GPS, what = "predicted")
PTT_ssm_tracks <- grab(fit_PTT, what = "predicted")
GLS_ssm_tracks <- grab(fit_GLS, what = "predicted")

#combine all together
ssm_tracks <- rbind(get0("GPS_ssm_tracks"),
                    get0("PTT_ssm_tracks"),
                    get0("GLS_ssm_tracks"))

#join device ids with individual ids
individuals <- tracks %>% 
  select(id, individual_id, device_id) %>% 
  group_by(id) %>% 
  summarise(individual_id = first(individual_id),
            device_id = first(device_id))
ssm_tracks <- ssm_tracks %>% left_join(individuals)

#select important columns and standardise names
ssm_tracks <- ssm_tracks %>% 
  ungroup() %>%
  select(individual_id, device_id, date, lon, lat, x.se, y.se) %>%
  rename(datetime = date, lon_se_km = x.se, lat_se_km = y.se)

#export
saveRDS(ssm_tracks,
        file = paste0(species_code, "_state_space_modelled.RDS"))
