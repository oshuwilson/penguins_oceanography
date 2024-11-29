#automated script for any seabird tracking data

rm(list=ls())

{
  library(tidyverse)
  library(hms)
}


##----------------------------------------------##
## User Input Only Required for start of Step 1 ##
##----------------------------------------------##

# 1. Read in and format track dataframe

#set required variables not in the data
species_code <- "ADPE"
study_code <- "NPOLAR_KG"
data_contact <- "Andy Lowther"
contact_email <- "andrew.lowther@npolar.no"
track_device <- "PTT"

##----------------##
## User Input End ##
##----------------##

#check that study hasn't already been processed before beginning
all_meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
if(study_code %in% all_meta$dataset_identifier){
  rm(list=ls())
  stop("The metadata for this study has already been processed.")
}

#set working directory to read files from
setwd(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code, "/Data"))

#list all available files
files <- list.files(pattern="*.csv")

#read in tracks
tracks <- read.csv(files)

#format columns
tracks <- tracks %>%
  mutate(datetime = paste(date_gmt, time_gmt, sep = " "),
         individual_id = as.factor(bird_id),
         device_id = as.factor(track_id),
         deployment_site = paste(colony_name, site_name, sep = ", ")) %>%
  mutate(datetime = as_datetime(datetime, format = "%Y-%m-%d %H:%M:%S")) %>%
  rename(lat = latitude, 
         lon = longitude,
         age_class = age,
         device_type = device,
         deployment_decimal_latitude = lat_colony,
         deployment_decimal_longitude = lon_colony) %>%
  arrange(individual_id, device_id, datetime)


# 2. Format metadata

#create metadata from tracking data
meta <- tracks %>% 
  group_by(individual_id, device_id) %>%
  summarise(dataset_identifier = study_code, abbreviated_name = species_code, 
            data_contact = data_contact, contact_email = contact_email,
            file_name = files, keepornot = NA, 
            scientific_name = first(scientific_name), common_name = first(common_name),
            deployment_site = first(deployment_site), age_class = first(age_class),
            sex = first(sex), device_type = first(device_type),
            deployment_decimal_latitude = first(deployment_decimal_latitude),
            deployment_decimal_longitude = first(deployment_decimal_longitude)) %>%
  ungroup()

#extract deployment date and times
deployments <- tracks %>%
  group_by(individual_id, device_id) %>%
  arrange(individual_id, device_id, datetime) %>%
  summarise(deployment_datetime = first(datetime)) %>%
  ungroup() %>%
  mutate(deployment_year = year(deployment_datetime),
         deployment_month = month(deployment_datetime),
         deployment_day = day(deployment_datetime),
         deployment_time = hms::as_hms(deployment_datetime)) %>%
  mutate(deployment_time = as.character(deployment_time))

#join deployment times to metadata
meta <- meta %>%
  left_join(select(deployments, individual_id, device_id, deployment_year, 
                   deployment_month, deployment_day, deployment_time))

#join to existing metadata
all_meta <- all_meta %>% bind_rows(meta)

#export metadata
saveRDS(meta, "metadata.RDS")
saveRDS(all_meta, "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")


# 3. Standardise Tracking Data

#create a new column for argos quality
tracks <- tracks %>%
  mutate(lc = if_else(device_type == "PTT", argos_quality,
                      if_else(device_type == "GPS", "G", "GL")))

#select the six required columns
tracks <- tracks %>%
  rename(date = datetime) %>%
  select(device_id, individual_id, date, lon, lat, lc)

#remove no data values and NAs
no_data_vals <- c(0, "LAT", "LON") #common no_data_vals include 0, LAT, LON
tracks <- tracks %>% 
  filter(!lat %in% no_data_vals & !lon %in% no_data_vals) %>%
  filter(!is.na(lat) & !is.na(lon))

#filtering
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/functions/standardise.R")
tracks <- standardise(df = tracks, #tracking dataframe
                      device_type = track_device) #either PTT, GPS, or GLS

#check for missing data on leap years
#extract tracks from leap years
leap_tracks <- tracks %>%
  filter(leap_year(date) == TRUE)

#extract data from Feb 28th, 29th, and March 1st
leap_dates <- leap_tracks %>%
  filter(month(date) == 2 & day(date) %in% c(28, 29) | month(date) == 3 & day(date) == 1)

#subset to each date
feb28 <- filter(leap_tracks, yday(date) == 59)
feb29 <- filter(leap_tracks, yday(date) == 60)
mar01 <- filter(leap_tracks, yday(date) == 61)

#if there is no data on Feb 29th but data on Feb 28th and March 1st, stop and warn user
if(nrow(feb29) == 0 & nrow(feb28) > 0 & nrow(mar01) > 0){
  stop("There is no data on February 29th for leap years. Please check these tracks before proceeding")
}

#cleanup
rm(leap_tracks, leap_dates, feb28, feb29, mar01)

#export standardised data
saveRDS(tracks, paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code, "/standardised_tracks.RDS"))
