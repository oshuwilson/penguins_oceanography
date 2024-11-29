#standardisation

#cleanup
rm(list=ls())

{
  library(data.table)
  library(tidyverse)
  library(terra)
}

##-----------------------------------------##
## User Input Only Required for Step 1 + 2 ##
##-----------------------------------------##

# 1. Read in tracking data and relevant metadata

#set required variables
species_code <- "ADPE" #the 4-letter RAATD code
study_code <- "ADPE_CEBC_PG" #Dataset identifier
device_type <- "GPS" #either PTT, GLS, or GPS
zonediff <- 0 #timezone difference relative to UTC

#set working directory to draw raw files from
setwd(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code, "/Data"))

#list all available files
files <- list.files(pattern="*.csv") #change pattern depending on filetype

#are all tracks in one file?
onefile <- if_else(length(files) == 1, T, F)

#read tracks - if all IDs in one file
if(onefile == TRUE){
  tracks <- read.csv(files)
}

#read tracks - if IDs are in separate files
#make sure to change substring numbers
if(onefile == FALSE){
  tracks <- NULL
  for(i in files){
    fixes <- read.csv(i)
    fixes$individual_id <- substr(i, 1, 10) #assign id from filenames
    fixes <- fixes %>% mutate_if(is.integer, as.character)
    tracks <- bind_rows(tracks, fixes)
    rm(fixes)
  }
}

#format date column as POSIXct date-time object
tracks <- tracks %>%
  mutate(date = as_datetime(datetime))


# 2. Format data for prefiltering

#create lc column for GPS or GLS data
if(device_type == "GPS"){
  tracks$lc <- "G"
}

if(device_type == "GLS"){
  tracks$lc <- "GL"
}

#read in device_id from metadata
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% filter(dataset_identifier == study_code & abbreviated_name == species_code) %>%
  mutate(individual_id = as.character(individual_id))

#join together - ADAPT FOR MULTIPLE DEVICES PER INDIVIDUAL
tracks <- tracks %>% 
  mutate(individual_id = as.character(individual_id)) %>%
  rename() %>%
  left_join(select(meta, individual_id, device_id))

#rename and select the six required columns
tracks <- tracks %>%
  rename(lon = lon,
         lat = lat,
         lc = lc) %>%
  select(device_id, individual_id, date, lon, lat, lc)

#make individual_id and device_id factors
tracks <- tracks %>% 
  mutate(device_id = as.factor(device_id),
         individual_id = as.factor(individual_id))

#change timezone
tracks$date <- tracks$date - hours(zonediff)

#remove no data values and NAs
no_data_vals <- c(0, "LAT", "LON", NA) #common no_data_vals include 0, LAT, LON
tracks <- tracks %>% 
  filter(!lat %in% no_data_vals & !lon %in% no_data_vals) %>%
  filter(!is.na(lat) & !is.na(lon))


##----------------##
## User Input End ##
##----------------##

# 3. Filtering
source("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Code/functions/standardise.R")

tracks <- tracks %>% mutate(lat = as.numeric(lat), lon = as.numeric(lon)) 

tracks <- standardise(df = tracks, #tracking dataframe
                      device_type = device_type) #either PTT, GPS, or GLS



# 4. Checks

# 4.1 check that redeployed device_ids don't have temporal overlap

#identify redeployed devices
redeployments <- meta %>%
  group_by(device_id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  filter(n > 1) %>%
  pull(device_id)

#extract first and last datetime for each device deployment
for(i in redeployments){
  #identify individual ids corresponding to this device id
  individual_ids <- meta %>% 
    mutate(individual_id = as.character(individual_id),
           device_id = as.character(device_id)) %>%
    filter(device_id == i) %>%
    pull(individual_id)
  
  #extract first and last datetime for each individual id
  for(j in individual_ids){
    #extract track data for that individual
    track <- tracks %>%
      filter(individual_id == j & device_id == i)
    
    #extract first and last datetime for that individual
    first <- track %>% 
      slice_head(n = 1) %>%
      pull(date)
    last <- track %>%
      slice_tail(n = 1) %>%
      pull(date)
    
    #create dataframe with first and last datetime
    if(j == individual_ids[1]){
      deployment_dates <- data.frame(individual_id = j,
                                     device_id = i,
                                     first = first,
                                     last = last)
    } else {
      deployment_dates <- deployment_dates %>% 
        bind_rows(data.frame(individual_id = j,
                             device_id = i,
                             first = first,
                             last = last))
    }
    
    #cleanup
    rm(track, first, last)
  }
  
  #check for overlap in deployment dates
  for(k in 1:nrow(deployment_dates)){
    
    #extract date range for deployment
    deployment_range <- c(deployment_dates$first[k], deployment_dates$last[k])
    
    #check whether other start or end dates fall into this range
    overlap <- deployment_dates %>%
      filter(individual_id != deployment_dates$individual_id[k]) %>%
      mutate(overlap = ifelse(first >= deployment_range[1] & first <= deployment_range[2] | 
                                last >= deployment_range[1] & last <= deployment_range[2], 
                              1, 0))
    
    #if there is overlap in any track, stop and warn user
    if(sum(overlap$overlap) > 0){
      stop(paste("There is temporal overlap of device_id", i, "for individual", deployment_dates$individual_id[k], 
                  "and individual(s)", paste(overlap$individual_id[overlap$overlap == 1], collapse = ", "), 
                 "\n Please check and amend these tracks before proceeding"))
    }
    
    #cleanup
    rm(deployment_range, overlap)
  }
}

#cleanup
rm(redeployments, individual_ids, deployment_dates)


# 4.2 Check that data isn't missing from February 29th on leap years

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


# 5. Export standardised data for trimming

#export standardised data
saveRDS(tracks, paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code, "/standardised_tracks.RDS"))
