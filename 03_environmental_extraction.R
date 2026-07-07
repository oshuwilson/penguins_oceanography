#-------------------------------------------------------------------------------
# Extract environmental variables to tracks and background samples
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
}

# prioritise terra extract
extract <- terra::extract

#-----------------------------------------------------------------------------
# Data Preparation
#-----------------------------------------------------------------------------

#define species
this.species <- "KIPE"

#load in species/region/stage info for this species
srs <- read.csv("data/tracks/species_site_stage_v2.csv")
srs <- srs %>% 
  filter(species == this.species)

#list each region
regions <- unique(srs$site)

#loop over all regions
for(this.site in regions){
  
  #identify stages for this region
  stages <- srs %>% 
    filter(site == this.site) %>% 
    pull(stage)
  
  #loop over each stage
  for(this.stage in stages){
    
    #load in tracks and background data
    tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))
    back <- readRDS(paste0("output/background/", this.species, "/", this.site, " ", this.stage, " background.RDS"))
    
    # reproject tracks to epsg:4326
    tracks <- tracks %>%
      vect(geom = c("x", "y"),
           crs = "epsg:6932") %>%
      project("epsg:4326") %>%
      as.data.frame(geom = "XY") %>%
      mutate(region = this.site, pa = "presence")
    
    # add presence-absence column to background data
    back <- back %>%
      mutate(region = this.site, pa = "absence")
    
    #combine the two datasets together
    data <- bind_rows(tracks, back)
    
    #only keep important columns
    data <- data %>%
      select(individual_id, date, deployment_site,
             region, pa, state, x, y)
    
    #cleanup
    rm(tracks, back)
    
    
    #-----------------------------------------------------------------------------
    # Extract static variables
    #-----------------------------------------------------------------------------
    
    #depth
    depth <- rast("E:/Satellite_Data/static/depth/depth.nc")
    
    #create SpatVector for tracks and background
    data <- vect(data,
                 geom=c("x", "y"),
                 crs=crs(depth)) #this ensures crs are the same as rasters
    
    #extract
    data$depth <- extract(depth, data, ID=F)
    
    #remove rows where depth is NA - will be NA for every GLORYS variable
    data <- data %>% drop_na(depth)
    
    
    #-----------------------------------------------------------------------------
    # Extract dynamic variables
    #-----------------------------------------------------------------------------
    
    #load in dynamic_extract functions
    source("code/functions/extraction_functions.R")
    
    #sea ice concentration
    data <- dynamic_extract("sic", data)
    data$sic[is.na(data$sic) & year(data$date) > 1992] <- 0 #SIC values of 0 print as NA in GLORYS
    print("sic")
    
    #current velocity
    data <- dynamic_extract("uo", data)
    data <- dynamic_extract("vo", data)
    data$curr <- sqrt((data$uo^2) + (data$vo^2))
    print("curr")
    
    #eddies
    data <- dynamic_extract("eddies", data)
    print("eddies")
    
    
    #-----------------------------------------------------------------------------
    # Export data 
    #-----------------------------------------------------------------------------
    
    #convert to data frame
    data <- as.data.frame(data, geom = "XY")
    
    #export
    saveRDS(data, file = paste0("output/extractions/", this.species, "/", this.site, " ", this.stage, " extracted.RDS"))
    
    #print completion
    print(paste0("Extraction complete for ", this.species, " ", this.site, " ", this.stage))
    
  }
}

