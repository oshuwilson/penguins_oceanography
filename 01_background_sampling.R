#-------------------------
#create background samples 
#-------------------------

#cleanup
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(terra)
  library(tidyverse)
}

#set species
this.species <- "CHPE"

#read in oceans and coast files for masking
oceans <- readRDS("data/oceans_vect.RDS")
coast <- readRDS("data/coast_vect.RDS")

#read in species/region/stage info 
srs <- read.csv("data/tracks/species_site_stage.csv")

#filter to this species
srs <- srs %>% 
  filter(species == this.species)

#list each region
regions <- unique(srs$site)

#loop over all regions
for(i in 1:length(regions)){
  
  #define region
  this.region <- regions[i]
  
  #list each stage name
  stages <- srs %>% 
    filter(site == this.region) %>% 
    pull(stage)
  
  #for each stage
  for(j in stages){

    #read in tracks
    tracks <- readRDS(paste0("output/tracks/", this.species, "/", this.region, " ", j, " tracks.RDS"))
    
    #convert to terra
    tracks_terra <- tracks %>%
      vect(geom = c("lon", "lat"),
           crs = "epsg:4326")
    
    #if lon spans -180/180, project to stereographic view
    if(min(tracks$lon) < -175 & max(tracks$lon) > 175){
      tracks_terra <- project(tracks_terra, "epsg:6932")
    }
    
    #visualise
    plot(tracks_terra, pch =".")
    
    #create minimum convex hull
    mch <- convHull(tracks_terra)
    plot(mch, add=T)
    
    #buffer to prevent self intersection error
    mch <- buffer(mch, 0)
    
    #if in polar projection, use coastfile for masking
    if(min(tracks$lon) < -175 & max(tracks$lon) > 175){
      mask <- erase(mch, coast)
    } else {
      mask <- terra::intersect(mch, oceans)
    }
    plot(mask)

    #create background points
    back <- spatSample(mask, size = nrow(tracks))
    plot(back, pch = ".", add = T)
    
    #convert to dataframe
    if(min(tracks$lon) < -175 & max(tracks$lon) > 175){
      back <- project(back, "epsg:4326")
    } 
    back <- as.data.frame(back, geom = "XY")
    back <- back %>% select(x, y)
    
    #link with a tracking ID and date
    back <- back %>% 
      mutate(individual_id = as.character(tracks$individual_id),
             date = tracks$date,
             region = this.region)
    
    #export 
    saveRDS(back, file = paste0("output/background/", this.species, "/", this.region, "_", j, "_background.RDS"))
    
    #print completion
    print(paste0(this.region, " ", j, " completed"))
  }
}
