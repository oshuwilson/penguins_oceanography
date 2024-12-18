#-------------------------------
#extract environmental variables
#-------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
}


# 1. Setup

#define species
this.species <- "ADPE"

#load in species/region/stage info for this species
srs <- read.csv("data/tracks/species_site_stage.csv")
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
    tracks <- readRDS(paste0("output/tracks/", this.species, "/", this.site, " ", this.stage, " tracks.RDS"))
    back <- readRDS(paste0("output/background/", this.species, "/", this.site, "_", this.stage, "_background.RDS"))
    
    #combine the two datasets together
    tracks <- tracks %>%
      rename(x = lon, y = lat) %>%
      mutate(region = this.site, 
             pa = "presence")
    back <- back %>%
      mutate(pa = "absence")
    data <- bind_rows(tracks, back)
    
    #cleanup
    rm(tracks, back)
    
    
    # 2. Extract environmental variables
    
    # 2.1 Static Variables
    
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
    
    #slope
    slope <- rast("E:/Satellite_Data/static/slope/slope.nc")
    data$slope <- extract(slope, data, ID=F)
    
    #dshelf
    dshelf <- rast("E:/Satellite_Data/static/dshelf/dshelf_resampled.nc")
    data$dshelf <- extract(dshelf, data, ID=F)
    
    #cleanup static
    rm(depth, slope, dshelf)
    
    
    # 2.2 Dynamic Variables 
    
    #load in dynamic_extract functions
    source("code/functions/extraction_functions.R")
    
    #sst 
    data <- dynamic_extract("sst", data)
    print("sst")
    
    #mld
    data <- dynamic_extract("mld", data)
    print("mld")
    
    #sal
    data <- dynamic_extract("sal", data)
    print("sal")
    
    #ssh
    data <- dynamic_extract("ssh", data)
    print("ssh")
    
    #sic
    data <- dynamic_extract("sic", data)
    data$sic[is.na(data$sic) & year(data$date) > 1992] <- 0 #SIC values of 0 print as NA in GLORYS
    print("sic")
    
    #curr
    data <- dynamic_extract("uo", data)
    data <- dynamic_extract("vo", data)
    data$curr <- sqrt((data$uo^2) + (data$vo^2))
    print("curr")
    
    
    # 2.3 Oceanographic/Cryospheric Variables
    
    #front_freq
    data <- dynamic_extract("front_freq", data)
    print("front_freq")
    
    #eddies
    data <- dynamic_extract("eddies", data)
    print("eddies")
    
    #dist2ice
    data <- dynamic_extract("dist2ice", data)
    print("dist2ice")
    
    #leads
    leads <- rast("E:/Satellite_Data/static/leads/leads_resampled.nc")
    data$leads <- extract(leads, data, ID=F)
    data <- data %>%
      mutate(leads = ifelse(month(date) %in% c(4, 5, 6, 7, 8, 9, 10), leads, NA))
    rm(leads)
    print("leads")
    
    #polynyas
    warp <- project(data, "epsg:3412")
    warp <- dynamic_extract("polynyas", warp)
    data <- project(warp, crs(data))
    print("polynyas")
    
    
    # 3. Format for export
    
    #convert to data frame
    data <- as.data.frame(data, geom = "XY")
    
    #export
    saveRDS(data, file = paste0("output/extractions/", this.species, "/", this.site, "_", this.stage, "_extracted.RDS"))
    
    #print completion
    print(paste0("Extraction complete for ", this.species, " ", this.site, " ", this.stage))
    
  }
}

