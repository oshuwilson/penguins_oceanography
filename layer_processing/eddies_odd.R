rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/")

{
  library(terra)
  library(tidyterra)
  library(ncdf4)
  library(lubridate)
  library(dplyr)
}

#------------ --------#
## Set Up Daily Loop ##
#------- -------------#

#read in cyclonic and anticyclonic eddies
all_cyclones <- readRDS("eddies/all_cyclones.RDS")
all_anticyclones <- readRDS("eddies/all_anticyclones.RDS")

#loop over each year
for(z in 2001){
  
  #subset to only cyclones and anticyclones from that year
  yearly_cyclones <- all_cyclones %>%
    filter(year(time) == z)
  yearly_anticyclones <- all_anticyclones %>%
    filter(year(time) == z)
  
  #set max number of days in year for loop
  maxdays <- if_else(leap_year(z), 366, 365)
  
  #loop over each day
  for(j in 1:maxdays){
    
    #filter cyclones and anticylones to this day
    cyclones <- yearly_cyclones %>%
      filter(yday(time) == j)
    anticyclones <- yearly_anticyclones %>%
      filter(yday(time) == j)
    
    #create ID column for later
    cyclones$id <- 1:length(cyclones$radius)
    anticyclones$id <- 1:length(anticyclones$radius)
    
    #form terra object for centers
    cyclones <- vect(cyclones,
                     geom = c("lon", "lat"),
                     crs="epsg:4326")
    cyclones <- project(cyclones, "EPSG:6932")
    
    anticyclones <- vect(anticyclones,
                         geom = c("lon", "lat"),
                         crs="epsg:4326")
    anticyclones <- project(anticyclones, "EPSG:6932")
    
    
    #---------------------------------#
    ## Create Cyclone Distance Layer ##
    #---------------------------------#
    
    #create 2 x Radius polygons
    double <- 2 * cyclones$radius
    doubles <- buffer(cyclones, width=double)
    
    #import depth layer to create raster
    depth <- rast("depth/depth_stereographic_2.RDS")
    e <- ext(doubles)
    depth <- crop(depth, e)
    
    #null for loop
    all_rel_dists <- NULL
    
    #create distance layer for each eddy
    for(i in 1:length(cyclones)){
      
      #remove previous loop variables
      rm(this_center, this_double, this_dist, this_radius, this_rel_dist)
      
      #compute distance to eddy center
      this_center <- cyclones[i]
      this_center <- rasterize(this_center, depth)
      this_center <- distance(this_center)
      this_center <- this_center %>% rename(eddy_dist = last)
      
      #mask out double-radius zone
      this_double <- doubles[i]
      this_dist <- mask(this_center, this_double)
      
      #calculate relative distance
      this_radius <- rasterize(this_double, depth, field="radius")
      this_rel_dist <- this_dist/this_radius
      
      #bind with other eddies
      all_rel_dists[[i]] <- this_rel_dist
      
      
    }
    
    #join each relative distance into one raster, with the minimum distance taking priority
    complete_dists <- rast(all_rel_dists)
    complete_dists <- min(complete_dists, na.rm=T) + 1
    cyclone_dists <- complete_dists
    
    
    #-------------------------------------#
    ## Create Anticyclone Distance Layer ##
    #-------------------------------------#

    #create 2 x Radius polygons
    double <- 2 * anticyclones$radius
    doubles <- buffer(anticyclones, width=double)
    
    #null for loop
    all_rel_dists <- NULL
    
    #create distance layer for each eddy
    for(i in 1:length(anticyclones)){
      this_center <- anticyclones[i]
      this_center <- rasterize(this_center, depth)
      this_center <- distance(this_center)
      this_center <- this_center %>% rename(eddy_dist = last)
      
      #mask out double-radius zone
      this_double <- doubles[i]
      this_dist <- mask(this_center, this_double)
      
      #calculate relative distance
      this_radius <- rasterize(this_double, depth, field="radius")
      this_rel_dist <- this_dist/this_radius
      
      all_rel_dists[[i]] <- this_rel_dist
    }
    
    #join each relative distance into one raster, with the minimum distance taking priority
    complete_dists <- rast(all_rel_dists)
    complete_dists <- min(complete_dists, na.rm=T) + 1
    anticyclone_dists <- complete_dists
    
    
    #--------------------#
    ## Join Eddy Layers ##
    #--------------------#
    
    #join layers and choose the minimum relative distance where they overlap
    eddy_dists <- rast(list(cyclone_dists, anticyclone_dists))
    eddy_dists <- min(eddy_dists, na.rm=T)
    
    #mask out cyclone and anticyclone layers so that only matching values are kept
    m <- c(-Inf, 0, 1,
           0, Inf, NA)
    m <- matrix(m, ncol=3, byrow=T)
    
    #cyclones
    cyclone_mask <- cyclone_dists - eddy_dists
    cyclone_mask <- classify(cyclone_mask, m)
    cyclone_mask <- as.polygons(cyclone_mask)
    cyclone_dists_2 <- mask(cyclone_dists, cyclone_mask)
    cyclone_dists_2 <- cyclone_dists_2 * -1 + 4
    
    #anticyclones
    anticyclone_mask <- anticyclone_dists - eddy_dists
    anticyclone_mask <- classify(anticyclone_mask, m)
    anticyclone_mask <- as.polygons(anticyclone_mask)
    anticyclone_dists_2 <- mask(anticyclone_dists, anticyclone_mask)
    anticyclone_dists_2 <- anticyclone_dists_2 - 4
    
    #rejoin together
    eddy_dists_2 <- rast(list(cyclone_dists_2, anticyclone_dists_2))
    eddy_dists_2 <- sum(eddy_dists_2, na.rm=T)
    
    #final raster for that day
    eddy_dists_2 <- subst(eddy_dists_2, NA, 0)
    
    #append time info
    time(eddy_dists_2) <- first(cyclones$time)
    
    #reproject to GLORYS resolution using SSH raster
    ssh <- rast("ssh/ssh_1993.nc")
    eddy_dists_2 <- project(eddy_dists_2, crs(ssh))
    
    #resample to SSH resolution
    eddy_dists_2 <- resample(eddy_dists_2, ssh, method = "bilinear")
    
    #crop to below 40 degrees south
    eddy_dists_2 <- crop(eddy_dists_2, ext(-180, 180, -80, -40))
    
    #join to raster for every day
      if(j == 1){
        yearly_eddies <- eddy_dists_2
      }
      
      if(j != 1){
        yearly_eddies <- c(yearly_eddies, eddy_dists_2)
      }
    
    #print that day has completed
    print(paste0(j, "/", maxdays, " completed"))
    
  }
  
  #export yearly raster
  writeCDF(yearly_eddies,
           filename = paste0("eddies/eddies_", z, ".nc"),
           varname = "eddies",
           longname = "relative distance to the cyclonic/anticyclonic eddy center")
  
  #print that year has completed
  print(paste0(z, " completed"))
  
}