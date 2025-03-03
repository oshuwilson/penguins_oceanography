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
this.species <- "MAPE"
this.study <- "NMU_Marion_2018"

#read in ssm_tracks and metadata for this species and study
tracks <- readRDS(paste0("data/tracks/", this.species, "/all_tracks.RDS"))
meta <- readRDS("data/metadata.RDS")
meta <- meta %>% 
  filter(abbreviated_name == this.species & 
           dataset_identifier == this.study)

#filter tracks to this study
tracks <- tracks %>% 
  filter(individual_id %in% meta$individual_id &
           device_id %in% meta$device_id)

#get each colony name
colony_names <- meta %>%
  select(deployment_site) %>%
  distinct() %>%
  pull()

#remove Marion Island colony name
colony_names <- colony_names[-which(colony_names == "Marion Island, Prince Edward Islands")]

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
  print(plot(trax, pch = "."))
  
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

#new early chick-rearing IDs
ecr <- c("1_2", "10_2", "10_3", "10_4", "10_5", "10_6", "10_7", "12_2", 
         "13_2", "13_3", "13_4", "13_5", "15_2", "2_2", "2_3", "2_4", "2_5", "2_6", "2_7",
         "34_2", "43_1", "44_1", "5_2", "5_3", "5_4", "5_5", 
         "8_2", "8_3", "8_4", "8_5", "8_6", "8_7", "8_8", "8_9", "9_2", "9_3",
         "17_2", "17_3", "18_2", "18_3", "18_4", "18_5", "18_6", "18_7", "18_8", "18_9", "18_10",
         "18_11", "19_2", "19_3", "19_4", "20_2", "21_3", "21_4", "21_5", "22_2", "22_3", "22_4",
         "22_5", "22_6", "22_7", "23_3", "23_4", "23_5", "23_6", "23_7", "23_8", "23_9", "23_10",
         "23_11", "24_2", "24_3", "24_4", "24_5", "24_6", "24_7", "24_8", "24_9", "24_10", "24_11",
         "24_12", "24_13", "24_14", "26_2", "26_3", "26_4", "27_2", "27_3", "27_4", "27_5", "27_6",
         "27_7", "27_8", "27_9", "27_10", "27_11", "27_12", "41_2", "41_3", "41_4", "41_5", "41_6",
         "42_2", "42_3")

# new incubation IDs
inc <- c("1_1", "10_1", "11_1", "12_1", "13_1", "14_1", "15_1", "2_1", "3_1",
         "30_1", "31_1", "32_1", "33_1", "34_1", "35_1", "4_1", "5_1", "6_1", "7_1",
         "8_1", "9_1", "16_1", "24_1", "25_1", "26_1", "36_1", "42_1")

# new late chick-rearing IDs
lcr <- c("17_4", "21_6", "27_16", "27_17")

#reassign stage names
all_tracks <- all_tracks %>%
  mutate(stage = case_when(
    trip %in% ecr ~ "early chick-rearing",
    trip %in% inc ~ "incubation",
    trip %in% lcr ~ "late chick-rearing",
    TRUE ~ stage
  ))

stage_by_id <- all_tracks %>% 
  group_by(individual_id, device_id, stage) %>%
  summarise(start = first(date), end = last(date)) 


# 6. Create stage dates for each id

#join to all stages
stages <- readRDS(paste0("data/stages/", this.species, "_stages.RDS"))
stages <- stages %>%
  filter(!individual_id %in% stage_by_id$individual_id) %>%
  bind_rows(stage_by_id) %>%
  select(-device_id)


#save stage dates
saveRDS(stages, paste0("data/stages/", this.species, "_stages.RDS"))
