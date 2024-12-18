#change deployment lat/lon code

library(CCAMLRGIS)
library(terra)

coast <- vect(load_Coastline())

rm(list=setdiff(ls(), "coast"))

# first part of metadata processing - wait for prompts
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/metadata_pt1.R")

#format initial columns - datetime, individual_id, and lat/lon
tracks <- tracks %>% select(Datetime, Latitude, Longitude, individual_id)
tracks$datetime <- as_datetime(tracks$Datetime, format = "%d/%m/%Y %H:%M")
tracks$individual_id <- as.factor(tracks$individual_id)
tracks <- tracks %>% rename(lat = Latitude, lon = Longitude)

# second part of metadata processing
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/metadata_pt2.R")

#format data from source metadata

#relevant columns often found in metadata: 
# individual_id - REQUIRED
# device_id - REQUIRED
# device_type
# sex
# age_class
# deployment_site
# deployment_date - must be date or POSIT format
# deployment_time
# deployment_latitude
# deployment_longitude

#keep and rename relevant columns
ext_meta <- ext_meta %>% 
  rename(device_id = GPS) %>%
  mutate(individual_id = as.factor(GPS_FILE_NAME)) %>%
  select(device_id, individual_id)

#if device_id not available, define as study code followed by numbers
if(!exists("device_id", where = ext_meta)){
  ext_meta <- ext_meta %>% mutate(device_id = paste(individual_id, 1:nrow(ext_meta), sep = "_"))
}

#append device ID to tracks - edit for multiple devices per individual
tracks <- tracks %>%
  left_join(select(ext_meta, individual_id, device_id))

# third part of metadata processing
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/metadata_pt3.R")

#rename deployment_site in ext_meta to include beach/bay and island/region
if(exists("deployment_site", where = ext_meta)){
  ext_meta <- ext_meta %>%
    mutate(deployment_site = gsub("South Orkneys, Signy Island", "Signy Island, South Orkney", deployment_site))
}

# fourth part of metadata processing - if errors reassign meta as meta3_banked
#CHANGE DEPLOYMENT LAT LON SCRIPT
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/metadata_pt4.R")

#check for NA values in inappropriate meta columns
na_test <- meta %>% 
  filter(is.na(dataset_identifier) | is.na(file_name) | is.na(individual_id) | is.na(scientific_name) | 
           is.na(common_name) | is.na(abbreviated_name) |  is.na(device_id) | is.na(device_type) |
           is.na(deployment_site) | is.na(deployment_year) | is.na(deployment_month) | is.na(deployment_day) | 
           is.na(deployment_time) | is.na(deployment_decimal_latitude) | is.na(deployment_decimal_longitude) |
           is.na(data_contact) | is.na(contact_email))

if(nrow(na_test) > 0){
  stop("NA values found in inappropriate columns - check na_test df for NA locations")
}

#export
saveRDS(all_meta,
        file = "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
saveRDS(meta, "metadata.RDS")
