#-------------------------------------------------------
#create background samples for each post-breeding colony 
#-------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(geosphere)
}

# read in species, site, and stage info to loop over populations
srs <- read.csv("data/tracks/species_site_stage_v2.csv")

# keep stages that aren't central-place-foraging
srs <- srs %>% 
  filter(stage %in% c("post-breeding", "pre-moult", "post-moult", "fledglings") |
           species == "KIPE" & stage == "late chick-rearing") # KIPE late chick-rearing is free-roaming

# read in oceans and coast files for masking
oceans <- readRDS("data/oceans_vect.RDS")
coast <- readRDS("data/coast_vect.RDS")

# sites of interest
srs <- srs %>% 
  filter(site %in% c("Auster Rookery", "Pointe Geologie", "Taylor Glacier", "Admiralty Bay, South Shetland"))

# loop over colonies 
for(j in 1:nrow(srs)){
  
  # colony and breeding stage values
  this.species <- srs$species[j]
  this.site <- srs$site[j]
  this.stage <- srs$stage[j]
  
  # 1. Create Polygon to Sample Within
  
  # load in hmm checked tracks
  tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))
  
  # convert to terra
  tracks_terra <- tracks %>%
    vect(geom = c("x", "y"),
         crs = "epsg:6932") %>%
    project("epsg:4326")
  
  # if lon spans -180/180, project to stereographic view
  if(ext(tracks_terra)[1] < -175 & ext(tracks_terra)[2] > 175){
    tracks_terra <- project(tracks_terra, "epsg:6932")
  }
  
  # visualise
  plot(tracks_terra, pch =".")
  
  # create minimum convex hull
  mch <- convHull(tracks_terra)
  plot(mch, add=T)
  
  # buffer to prevent self intersection error
  mch <- buffer(mch, 0)
  
  # if in polar projection, use coastfile for masking
  if(ext(tracks_terra)[1] < -175 & ext(tracks_terra)[2] > 175){
    mask <- erase(mch, coast)
  } else {
    mask <- terra::intersect(mch, oceans)
  }
  plot(mask)
  
  
  # 2. Generate Background Samples
  
  # define number of samples as 20,000 if nrow(tracks) < 20,000
  if(nrow(tracks) < 20000){
    n <- 20000
  } else {
    n <- nrow(tracks)
  }
  
  # create background points
  back <- spatSample(mask, size = n)
  plot(back, pch = ".", add = T)
  
  # assign dates and individual IDs proportionately from tracks
  back <- back %>%
    mutate(date = sample(tracks$date, n(), replace = T),
           individual_id = sample(tracks$individual_id, n(), replace = T))
  
  # reproject if in polar projection
  if(ext(tracks_terra)[1] < -175 & ext(tracks_terra)[2] > 175){
    back <- project(back, "epsg:4326")
  } 
  
  # 3. Extract Covariate Values
  
  # load dynamic extract function
  source("code/functions/extraction_functions.R")
  
  # extract eddies 
  back <- dynamic_extract("eddies", back)
  
  # extract uo and vo to calculate current speed
  back <- dynamic_extract("uo", back)
  back <- dynamic_extract("vo", back)
  back <- back %>%
    mutate(curr = sqrt(uo^2 + vo^2))
  
  # extract depth
  depth <- rast("E:/Satellite_Data/static/depth/depth.nc")
  back$depth <- extract(depth, back, ID=F)
  
  # for certain species, extract sea ice concentration
  if(this.species %in% c("ADPE", "CHPE", "EMPE")){
    back <- dynamic_extract("sic", back)
  }
  
  
  # 4. Export
  
  #convert to dataframe
  back <- as.data.frame(back, geom = "XY")
  
  # save
  saveRDS(back, paste0("output/background/", this.species, "/", this.site, " ", this.stage, " background.rds"))
  
  # print completion
  print(paste0(this.species, " ", this.site, " ", this.stage, " complete"))
}
