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
coast <- vect(load_Coastline())

#define species
this.species <- "ADPE"

#read in ssm_tracks and metadata for this species
tracks <- readRDS(paste0("data/ssm_tracks/", this.species, "_ssm_qc.RDS"))
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% 
  filter(abbreviated_name == this.species)


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

#convert tracks to terra stereographic view
trax <- vect(tracks, geom = c("lon", "lat"), crs = "epsg:4326")
trax <- project(trax, "epsg:6932")
plot(trax, pch = ".")

#split up trips
tracks <- trip_split(trax, meta, buff.dist = 5000)

#crop coastline to track extent
crop_coast <- crop(coast, ext(trax))


# 3. Assign stage dates to each trip

#read in stage dates
stage_dates <- readRDS(paste0("~/OneDrive - University of Southampton/Documents/Chapter 02/code/stage_splitting/config/", this.species, "config.RDS"))

#extract first point of each trip
trip_starts <- tracks %>%
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
tracks <- tracks %>%
  left_join(select(trip_starts, individual_id, device_id, trip, stage))

#cleanup
rm(trip_starts, stage_dates, this.stage, new_year, i)


# 4. Visualize tracks for each stage

#convert trip to factor
tracks$trip <- as.factor(tracks$trip)

#define all stages present
stages <- tracks %>%
  mutate(stage = as.factor(stage)) %>%
  pull(stage)
stages <- levels(stages)

#define all ids present
ids <- levels(tracks$id)

#setup pdf export
pdf(file = paste0("code/stage_splitting/checks/", this.species, "_checks.pdf"), width = 8, height = 10, pointsize = 16)

#plot tracks for each id and stage and highlight trips
for(i in ids){
  
  #subset tracks for this id
  trax <- tracks %>% filter(id == i)
  
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
    all_stage <- tracks %>% 
      filter(id != i & stage == j) %>%
      vect(geom = c("x", "y"), crs = "epsg:6932")
    
    #plot together
    p1 <- ggplot() + geom_spatvector(data = crop_coast, fill = "white") +
      geom_spatvector(data = all_stage, size = 0.5, color = "grey") +
      geom_spatvector(data = stage_trax, size = 2, shape = 17, aes(color = trip)) + 
      theme_bw() + scale_color_viridis_d() +
      labs(title = paste(i, j))
    p1
    print(p1)
  }
}

#finish pdf export
dev.off()


# 5. Create stage dates for each id (NEEDS MANUAL ADJUSTMENT STEP)

#format final stage dates
stage_by_id <- tracks %>% 
  group_by(individual_id, device_id, stage) %>%
  summarise(start = first(date), end = last(date))

#save initial stage dates
saveRDS(stage_by_id, paste0("code/stage_splitting/stage_dates/", this.species, "_stages.RDS"))

#save tracks with trips attached
saveRDS(tracks, paste0("code/stage_splitting/stage_dates/", this.species, "_tracks_with_stage_trips.RDS"))
