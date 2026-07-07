#-------------------------------------------------------------------------------
# Create background samples for each post-breeding case study 
#-------------------------------------------------------------------------------

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

# loop over colonies 
for(j in 1:nrow(srs)){
  
  # colony and breeding stage values
  this.species <- srs$species[j]
  this.site <- srs$site[j]
  this.stage <- srs$stage[j]
  
  #-----------------------------------------------------------------------------
  # Create minimum convex polygon
  #-----------------------------------------------------------------------------
  
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
  
  
  #-----------------------------------------------------------------------------
  # Generate background samples
  #-----------------------------------------------------------------------------
  
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
  
  
  #-----------------------------------------------------------------------------
  # Export outputs
  #-----------------------------------------------------------------------------

  #convert to dataframe
  back <- as.data.frame(back, geom = "XY")
  
  # save
  saveRDS(back, paste0("output/background/", this.species, "/", this.site, " ", this.stage, " background.rds"))
  
  # print completion
  print(paste0(this.species, " ", this.site, " ", this.stage, " complete"))
}
