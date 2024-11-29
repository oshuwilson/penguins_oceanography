#formatting RAATD metadata - run from source (Ctrl+Shift+S) NEED DIFFERENT META NAME FOR EACH STEP
{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(CCAMLRGIS)
  library(hms)
}

#coastfile for visualisations
coast <- vect(load_Coastline())
rm(list=setdiff(ls(), "coast"))


# 1. Read in and format track dataframe

#set species, study and colony/region code
species_code <- "ADPE"
study_code <- "ADPE_CEBC_PG"

#check that study hasn't already been processed before beginning
all_meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
if(study_code %in% all_meta$dataset_identifier){
  stop("The metadata for this study has already been processed.")
}

#set working directory to read files from
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
if(onefile == FALSE){
  first_num <- 
  last_num <- 
  tracks <- NULL
  for(i in files){
    fixes <- read.csv(i)
    fixes$individual_id <- substr(i, first_num, last_num) 
    tracks <- bind_rows(tracks, fixes)
    rm(fixes)
  }
}

#format initial columns - datetime, individual_id, and lat/lon
tracks$datetime <- as_datetime(tracks$datetime)
tracks$individual_id <- as.factor(tracks$individual_id)
tracks <- tracks %>% rename(lat = lat, lon = lon)

#change timezone to UTC
zonediff <- 0
tracks$datetime <- tracks$datetime - hours(zonediff)

#cleanup
rm(i, onefile, zonediff)


# 2. Read in Existing Metadata and keep relevant columns
#read in existing metadata
ext_meta <- read.csv(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code, "/metadata.csv"))

#relevant columns often found in metadata: 
# individual_id 
# device_id 
# device_type
# sex
# age_class
# deployment_site
# deployment_date
# deployment_time
# deployment_latitude
# deployment_longitude

#keep and rename relevant columns
ext_meta <- ext_meta %>% 
  select(New.GPS.name, GPS_no, sex, latitude_degS, longitude_degN) %>%
  rename(individual_id = New.GPS.name,
         device_id = GPS_no,
         deployment_latitude = latitude_degS,
         deployment_longitude = longitude_degN) 


#append device_id to tracks
tracks <- tracks %>% mutate(device_id = individual_id)


# 3. Create dataframe for individual study

### group by individual_id and device_id, and append study_code and species_code ###
meta <- tracks %>% group_by(individual_id, device_id) %>%
  summarise(dataset_identifier = study_code, abbreviated_name = species_code) %>%
  ungroup()

### scientific_name, common_name ### 
#import list of names
names <- read.csv("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/species_codes.csv")

#filter to this species
names <- filter(names, abbreviated_name == species_code)

#assign common and scientific name from this
meta$scientific_name <- names$scientific_name
meta$common_name <- names$common_name

#cleanup
rm(names)


### sex, age_class ###
#sex
#meta$sex <- "unknown"

meta <- meta %>% 
  left_join(select(ext_meta, individual_id, sex)) %>%
  distinct(individual_id, .keep_all = T)
meta <- meta %>%
  mutate(sex = case_match(sex,
         "M" ~ "male",
         "F" ~ "female", 
         "" ~ "unknown"))
meta <- meta %>% 
  mutate(sex = if_else(is.na(sex), "unknown", sex))

#check that sex is either male, female, or NA
poss_sexes <- c("male", "female", "unknown")
invalid <- meta %>%
  filter(!sex %in% poss_sexes)
if(nrow(invalid) > 0){
  print(invalid)
  stop("sex must be male, female or unknown. Please rename others.")
}

#cleanup
rm(invalid, poss_sexes)

#age_class - if in ext_meta
meta$age_class <- "adult"

#check that age_class is either adult, juvenile, or unknown
poss_ages <- c("adult", "juvenile", "unknown")
invalid <- meta %>%
  filter(!age_class %in% poss_ages)
if(nrow(invalid) > 0){
  print(invalid)
  stop("age_class must be adult, juvenile or unknown. Please rename others.")
}

#cleanup
rm(invalid, poss_ages)


### device_type ###

#device_type - if in ext_meta
meta$device_type <- "GPS"

#check that device_type has been saved as GPS, PTT, or GLS
poss_devices <- c("GPS", "PTT", "GLS")
invalid <- meta %>%
  filter(!device_type %in% poss_devices)
if(nrow(invalid) > 0){
  print(invalid)
  stop("device_type must be GPS, PTT, or GLS. Please rename others.")
}

#cleanup
rm(invalid, poss_devices)


### deployment_site ###
meta$deployment_site <- "Pointe Geologie"


### deployment date, time, and lat/lon ###
meta$deployment_decimal_latitude <- -66.66667
meta$deployment_decimal_longitude <- 140.01667

#if date, time, and/or lat/lons are not in ext_meta, use first point for each device_ID, but exercise caution
if(!exists("ext_meta")){
  ext_meta <- data.frame(0,0,0)
}

if(!exists("deployment_date", where = ext_meta) | !exists("deployment_lat", where = ext_meta)){
  tracks <- tracks %>% filter(lat != 0 & lon != 0 & lat != "LAT" & lon != "LON")
  tracks <- tracks %>% mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  
  deployments <- tracks %>% 
    group_by(individual_id, device_id) %>% 
    drop_na(all_of(c("lat", "lon"))) %>%
    summarise(start_time = first(datetime),
              start_lat = first(lat),
              start_lon = first(lon))
  
  #visualise to check
  deps <- vect(deployments, geom = c("start_lon", "start_lat"), crs = "EPSG:4326")
  deps <- project(deps, "EPSG:6932")
  e <- ext(deps) + c(10000, 10000, 10000, 10000)
  crop_coast <- crop(coast, e)
  plot(crop_coast)
  plot(deps, add = T)
  
  #extract deployment time and date
  deployments <- deployments %>% 
    mutate(deployment_year = year(start_time),
           deployment_month = month(start_time),
           deployment_day = day(start_time),
           deployment_time = as_hms(start_time))
  
  #format columns for metadata
  deployments <- deployments %>% 
    rename(deployment_decimal_latitude = start_lat,
           deployment_decimal_longitude = start_lon) %>%
    select(-start_time)
  
  #if happy, extract deployment information - only choose data missing from ext_meta (FIX THIS BIT)
  meta <- meta %>% 
    left_join(select(deployments, device_id, individual_id,
                     deployment_day, deployment_month, deployment_year, deployment_time),
              by = c("individual_id", "device_id"))
  
  #cleanup
  rm(deployments, deps, e, crop_coast)
}

### data contact and email ###
meta$data_contact <- "Akiko Kato"
meta$contact_email <- "akiko.k.r@gmail.com"


### split tracks into separate files and store filenames ###
#metadata filename column
meta$file_name <- files

### keepornot column for filtering ###
#NA for now - will be populated later
meta$keepornot <- NA


# 4. Join with existing metadata
#reorder columns
meta <- meta %>% select(dataset_identifier, file_name, individual_id, keepornot, scientific_name, 
                        common_name, abbreviated_name, sex, age_class, device_id, device_type, 
                        deployment_site, deployment_year, deployment_month, deployment_day,
                        deployment_time, deployment_decimal_longitude, deployment_decimal_latitude,
                        data_contact, contact_email)

#check whether there are duplicates before joining dataset
inds <- meta %>% 
  mutate(individual_id = as.character(individual_id)) %>%
  pull(individual_id)
devs <- meta %>% 
  mutate(device_id = as.character(device_id)) %>%
  pull(device_id)

dups <- all_meta %>% 
  filter(individual_id %in% inds & abbreviated_name == species_code |
           device_id %in% devs & abbreviated_name == species_code)

if(nrow(dups) > 0){
  print(dups)
  stop("Duplicate IDs found in metadata - check for temporal overlap in individuals or devices from dups dataset.")
}

all_meta <- all_meta %>%
  mutate(device_id = as.character(device_id))
meta <- meta %>%
  mutate(device_id = as.character(device_id))

meta <- meta %>%
  mutate(deployment_time = as.character(deployment_time))

#join together
all_meta <- all_meta %>% bind_rows(meta) %>%
  arrange(abbreviated_name, dataset_identifier, individual_id) 

#save updated metadata
saveRDS(all_meta,
        file = "~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
saveRDS(meta, file = "metadata.RDS")
