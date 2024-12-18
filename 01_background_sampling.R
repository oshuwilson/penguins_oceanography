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
this.species <- "ADPE"

#read in oceans file and metadata
oceans <- readRDS("data/oceans_vect.RDS")

#read in metadata and filter to this species
meta <- read.csv("~/OneDrive - University of Southampton/Documents/RAATD/RAATD_metadata.csv")
meta <- meta %>% 
  filter(abbreviated_name == this.species)

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
  
  #list stage track files
  files <- list.files(path = paste0("output/tracks/", this.species, "/"),
                      pattern = paste0(this.region),
                      full.names = T)
  
  #list each stage name
  stages <- srs %>% 
    filter(site == this.region) %>% 
    pull(stage)
  
  #for each stage
  for(j in 1:length(stages)){
    stage_file <- files[j]
    stage_name <- stages[j]
    
    #read in tracks
    tracks <- readRDS(stage_file)
    
    #convert to terra
    tracks_terra <- tracks %>%
      vect(geom = c("lon", "lat"),
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
    saveRDS(back, file = paste0("output/background/", this.species, "/", this.region, "_", stage_name, "_background.RDS"))
    
    #print completion
    print(paste0(this.region, " ", stage_name, " completed"))
  }
}
