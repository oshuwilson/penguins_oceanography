#create background samples for penguins

#RESAMPLE ALL TO A 2-HOUR STEP-LENGTH????
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(sf)
  library(terra)
  library(tidyverse)
}

#set species
this.species <- "ADPE"

#read in oceans file and metadata
oceans <- readRDS("data/oceans_vect.RDS")

#read in metadata and filter to this species
meta <- read.csv("~/OneDrive - University of Southampton/Documents/RAATD/RAATD_metadata.csv")
meta <- meta %>% 
  filter(abbreviated_name == this.species)

#create list of all regions
regions <- c("Pointe Geologie")

#loop over all regions
for(i in 1:length(regions)){
  
  #define region
  this.region <- regions[i]
  
  #list stage track files
  files <- list.files(path = paste0("data/tracks/", this.species, "/", this.region, "/"),
                      full.names = T)
  stages <- list.files(path = paste0("data/tracks/", this.species, "/", this.region, "/"))
  
  #for each stage
  for(j in 1:length(stages)){
    stage_file <- files[j]
    stage_name <- stages[j]
    stage_name <- tools::file_path_sans_ext(stage_name)
    
    #read in tracks
    tracks <- readRDS(stage_file)
    
    #convert to terra
    tracks_terra <- tracks %>%
      vect(geom = c("decimal_longitude", "decimal_latitude"),
           crs = "epsg:4326")
    
    #visualise
    plot(tracks_terra, pch =".")
    
    #create minimum convex hull
    mch <- convHull(tracks_terra)
    plot(mch, add=T)
    
    #buffer to prevent self intersection error
    mch <- buffer(mch, 0)
    
    #intersect mch with oceans to avoid sampling on land
    mask <- terra::intersect(mch, oceans)
    plot(mask)
    
    #create background points
    back <- spatSample(mask, size = nrow(tracks))
    plot(back, pch = ".", add = T)
    
    #convert to dataframe
    back <- as.data.frame(back, geom = "XY")
    back <- back %>% select(x, y)
    
    #link with a tracking ID and date
    back <- back %>% 
      mutate(individual_id = as.character(tracks$individual_id),
             date = tracks$date,
             region = this.region)
    
    #export 
    saveRDS(back, file = paste0("data/tracks/", this.species, "/", this.region, "/", stage_name, "_background.RDS"))
    
    #print completion
    print(paste0(this.region, " ", stage_name, " completed"))
  }
}
