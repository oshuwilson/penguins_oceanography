#split tracks into trips

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(CCAMLRGIS)
}

#read in coastline
coast <- readRDS("data/coast_vect.RDS")

#define species and study
this.species <- "EMPE"
this.study <- "Colbeck_2013"

#read in ssm_tracks and metadata for this species and study
tracks <- readRDS(paste0("data/tracks/", this.species, "/all_tracks.RDS"))
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% 
  filter(abbreviated_name == this.species & 
           dataset_identifier == this.study)

#change colony name for CHPE to avoid slash
meta <- meta %>%
  mutate(deployment_site = ifelse(deployment_site == "Signy Island/Gourlay, South Orkney Islands", 
                                  "Signy Island, South Orkney", deployment_site))

#filter tracks to this study
tracks <- tracks %>% 
  filter(individual_id %in% meta$individual_id &
           device_id %in% meta$device_id)

#get each colony name
colony_names <- meta %>%
  select(deployment_site) %>%
  distinct() %>%
  pull()

#read in stage dates
stage_dates <- readRDS(paste0("~/OneDrive - University of Southampton/Documents/Chapter 02/code/stage_splitting/config/", this.species, "config.RDS"))


# 1. Calculate the distance of locations from deployment sites by ID

#create combined ID
tracks <- tracks %>%
  mutate(id = as.factor(paste(individual_id, device_id, sep = "_")))
meta <- meta %>%
  mutate(id = as.factor(paste(individual_id, device_id, sep = "_")))


#null variable for loop
all_tracks <- NULL

#run over each ID
for(i in levels(tracks$id)){
  
  #tracks for this id
  trax <- tracks %>% filter(id == i)
  
  #deployment lat/lon for this id
  dep <- meta %>%
    filter(id == i) %>%
    select(deployment_decimal_longitude, deployment_decimal_latitude) %>%
    rename(lon = deployment_decimal_longitude, lat = deployment_decimal_latitude)
  
  #convert tracks and deployment to terra
  trax <- vect(trax, geom = c("lon", "lat"), crs = "epsg:4326")
  dep <- vect(dep, geom = c("lon", "lat"), crs = "epsg:4326")
  
  #calculate distance
  trax$distance <- distance(trax, dep)
  
  #convert back to dataframe
  trax <- as.data.frame(trax, geom="XY")
  trax <- trax %>%
    rename(lon = x, lat = y)
  
  #append to all tracks
  all_tracks <- bind_rows(all_tracks, trax)
  
}

#rename all_tracks and cleanup
tracks <- all_tracks
rm(all_tracks, dep, trax, i)


# 2. Split tracks into trips

#load trip splitting function
source("~/OneDrive - University of Southampton/Documents/Chapter 02/code/functions/trip_split.R")

#null for all tracks
all_tracks <- NULL

#loop over each colony
for(i in colony_names){
  
  colony_name <- i
  
  #split tracks to tracks in this colony
  colony_tracks <- tracks %>% 
    left_join(meta %>% select(individual_id, device_id, deployment_site)) %>%
    filter(deployment_site == i)
  
  #convert tracks to terra stereographic view
  trax <- vect(colony_tracks, geom = c("lon", "lat"), crs = "epsg:4326")
  trax <- project(trax, "epsg:6932")
  plot(trax, pch = ".")
  
  #split up trips
  colony_tracks <- trip_split(trax, meta, buff.dist = 10000)
  
  #crop coastline to track extent
  crop_coast <- crop(coast, ext(trax))
  
  
  # 3. Assign stage dates to each trip
  
  #extract first point of each trip
  trip_starts <- colony_tracks %>%
    group_by(id, trip) %>%
    slice(1) %>%
    ungroup() %>%
    select(individual_id, device_id, date, trip)
  
  #create yday column for trip starts
  trip_starts <- trip_starts %>%
    mutate(yday = yday(date))
  
  #create column for automatic stage assignment in trip starts
  trip_starts <- trip_starts %>%
    mutate(stage = NA)
  
  
  #join stage dates to trip starts by stage
  for(i in 1:nrow(stage_dates)){
    
    #define stage
    this.stage <- stage_dates[i,]
    
    #check if stage covers new year
    if(this.stage$end_day < this.stage$start_day){
      new_year <- TRUE
    } else {
      new_year <- FALSE
    }
    
    #assign stage to trip starts
    if(new_year == FALSE){
      trip_starts <- trip_starts %>%
        mutate(stage = ifelse(yday >= this.stage$start_day & yday <= this.stage$end_day, this.stage$stage, stage))
    } else {
      trip_starts <- trip_starts %>%
        mutate(stage = ifelse(yday >= this.stage$start_day | yday <= this.stage$end_day, this.stage$stage, stage))
    }
    
  }
  
  #join stage to tracks
  colony_tracks <- colony_tracks %>%
    left_join(select(trip_starts, individual_id, device_id, trip, stage))
  
  #cleanup
  rm(trip_starts, this.stage, new_year)
  
  
  # 4. Visualize tracks for each stage
  
  #convert trip to factor
  colony_tracks$trip <- as.factor(colony_tracks$trip)
  
  #define all stages present
  stages <- colony_tracks %>%
    mutate(stage = as.factor(stage)) %>%
    pull(stage)
  stages <- levels(stages)
  
  #define all ids present
  ids <- levels(colony_tracks$id)
  
  #setup pdf export
  pdf(file = paste0("code/stage_splitting/checks/", this.species, "_", colony_name, "_checks.pdf"), width = 8, height = 10, pointsize = 16)
  
  #plot tracks for each id and stage and highlight trips
  for(i in ids){
    
    #subset tracks for this id
    trax <- colony_tracks %>% filter(id == i)
    
    for(j in stages){
      
      #subset tracks for this stage
      stage_trax <- trax %>% filter(stage == j)
      
      #if empty skip
      if(nrow(stage_trax) == 0){
        next
      }
      
      #convert to terra
      stage_trax <- vect(stage_trax, geom = c("x", "y"), crs = "epsg:6932")
      
      #create terra object for all other tracks of this stage
      all_stage <- colony_tracks %>% 
        filter(id != i & stage == j) %>%
        vect(geom = c("x", "y"), crs = "epsg:6932")
      
      #plot together
      if(length(all_stage) > 0){
        p1 <- ggplot() + geom_spatvector(data = crop_coast, fill = "white") +
          geom_spatvector(data = all_stage, size = 0.5, color = "grey") +
          geom_spatvector(data = stage_trax, size = 2, shape = 17, aes(color = trip)) + 
          theme_bw() + scale_color_viridis_d() +
          labs(title = paste(i, j))
      } else{
        p1 <- ggplot() + geom_spatvector(data = crop_coast, fill = "white") +
          geom_spatvector(data = stage_trax, size = 2, shape = 17, aes(color = trip)) + 
          theme_bw() + scale_color_viridis_d() +
          labs(title = paste(i, j))
      }
      print(p1)
    }
  }
  
  #finish pdf export
  dev.off()
  
  #join colony tracks to all other tracks
  all_tracks <- bind_rows(all_tracks, colony_tracks)
  
}

#format final stage dates
stage_by_id <- all_tracks %>% 
  group_by(individual_id, device_id, stage) %>%
  summarise(start = first(date), end = last(date))


# 5. Reassign trips from visual checks or other info

#new chick-rearing IDs
chickrearing <- c("29330_17_12_13", "29330_27_12_13", "29602_20_12_13", 
                  "29606_18_12_13", "29626_18_12_13", "29629_18_12_13")

#reassign stage names
stage_by_id <- stage_by_id %>%
  mutate(stage = case_when(
    individual_id %in% chickrearing ~ "chick-rearing",
    TRUE ~ stage
  )) %>%
  ungroup() %>%
  group_by(individual_id, device_id, stage) %>%
  summarise(start = min(start), end = max(end))


# 6. Create stage dates for each id

#join to all stages
stages <- readRDS(paste0("data/stages/", this.species, "_stages.RDS"))
stages <- stages %>%
  filter(!individual_id %in% stage_by_id$individual_id &
           !device_id %in% stage_by_id$device_id) %>%
  bind_rows(stage_by_id)


#save stage dates
saveRDS(stages, paste0("data/stages/", this.species, "_stages.RDS"))
