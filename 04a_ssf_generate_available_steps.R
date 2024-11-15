#testing step-selection functions
#need to integrate device IDs?

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(amt)
  library(terra)
  library(tidyterra)
}

# 1. Setup

#test run with KIPE from Macquarie
this.species <- "KIPE"
this.site <- "Crozet"
this.stage <- "incubation"

#read in tracks
tracks <- readRDS(paste0("data/tracks/", this.species, "/", this.site, "/", this.stage, ".RDS"))
tracks <- tracks %>% 
  mutate(individual_id = as.character(individual_id)) %>%
  mutate(individual_id = as.factor(individual_id))

#read in metadata
meta <- read.csv("~/OneDrive - University of Southampton/Documents/RAATD/RAATD_metadata.csv")
meta <- meta %>% 
  filter(abbreviated_name == this.species,
         individual_id %in% levels(tracks$individual_id))

#project to planar projection
tracks <- vect(tracks, geom = c("decimal_longitude", "decimal_latitude"), crs = "epsg:4326")
tracks <- project(tracks, "epsg:6932")
plot(tracks, pch = ".")


# 2. Split tracks into trips 

#use function trip_split, use printed plot to check buffer distance is suitable
source("code/functions/trip_split.R")
tracks <- trip_split(tracks, meta, buff.dist = 10000)


# 3. Create available steps for each trip

#convert to amt format
amt_tracks <- make_track(tracks, x, y, date, 
                         crs = 4326, id = individual_id,
                         trip_id = trip)

#list of all individuals
trip_ids <- levels(as.factor(tracks$trip))

#null variable for loop
all_steps <- NULL

#loop over each trip
for(i in trip_ids){
  
  #filter to trip
  this.trip <- i
  trax <- amt_tracks %>% 
    filter(trip_id == this.trip)
  
  #resample to 6 hour steps - increases step length (will this cause issues in shorter foraging stages?)
  trax <- track_resample(trax, rate = hours(6), tolerance = minutes(30))
  
  #move on if resampled tracks are too small
  if(nrow(trax) < 4){
    next
  } 
  
  #convert track to steps
  steps <- steps_by_burst(trax)
  
  #create available steps
  steps <- steps %>% random_steps(n_control = 20)
  
  #append individual and trip ids from tracks
  steps$trip_id <- this.trip
  steps$individual_id <- first(trax$id)
  
  #create time from start of trip variable
  steps$time_since_start <- cumsum(as.numeric(steps$dt_))
  
  #calculate distance to colony
  source("code/functions/distance_to_colony.R")
  steps$dist_to_colony <- distance_to_colony(steps, meta)
  
  #append to null variable
  all_steps <- bind_rows(all_steps, steps)
  
}


# 4. Export

#save available and observed steps
saveRDS(all_steps, paste0("output/ssf/available steps/", this.species, "/", this.site, "_", this.stage, "_available_steps.RDS"))
