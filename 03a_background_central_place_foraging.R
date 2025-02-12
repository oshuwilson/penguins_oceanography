#---------------------------------------------------------------
# Background Sampling for Central Place Foraging Breeding Stages
#---------------------------------------------------------------

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

# remove stages that aren't central-place-foraging
srs <- srs %>% 
  filter(!stage %in% c("post-breeding", "pre-moult", "post-moult")) %>%
  filter(species != "KIPE" | stage != "late chick-rearing") #MAPE late chick-rearing is CPF

# loop over each breeding stage at each colony
for(j in 1:nrow(srs)){
  
  #cleanup 
  rm(list = setdiff(ls(), c("j", "srs")))
  
  # 1. Data Preparation
  
  # colony and breeding stage values
  this.species <- srs$species[j]
  this.site <- srs$site[j]
  this.stage <- srs$stage[j]
  
  # load in metadata for colony location
  meta <- readRDS("data/metadata.rds")
  
  # load in hmm checked tracks
  tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))
  
  # filter metadata to this colony
  meta <- meta %>% filter(abbreviated_name == this.species &
                            individual_id %in% tracks$individual_id)
  
  # convert to terra
  tracks <- tracks %>%
    vect(geom = c("x", "y"),
         crs = "epsg:6932") %>%
    project("epsg:4326")
  
  # load in coastline
  coast <- readRDS("data/land_vect.RDS")
  
  
  # 2. Create Buffer Around Colony
  
  # get colony location
  colony <- meta %>%
    slice(1) %>%
    select(lon = deployment_decimal_longitude, lat = deployment_decimal_latitude) %>%
    vect(geom = c("lon", "lat"), crs = "epsg:4326") 
  
  # calculate maximum distance to colony
  tracks$distances <- distance(tracks, colony)
  max_dist <- max(tracks$distances)
  
  # create buffer from colony of that distance
  buff <- buffer(colony, max_dist)
  
  
  # 3. Segment Buffer Based on Track Bearings from Colony
  
  # calculate bearings from colony
  trax <- tracks %>%
    filter(distances > 0.3 * max_dist) %>%
    crds() 
  col <- meta %>%
    slice(1) %>%
    select(lon = deployment_decimal_longitude, lat = deployment_decimal_latitude) %>%
    as.matrix()
  angles <- bearing(col, trax)
  
  # take min and max bearing, if cross the -180/180 line, add 360 to the negative values
  if(min(angles) < -178 & max(angles) > 178){
    angles <- ifelse(angles < 0, angles + 360, angles)
  }
  
  # calculate min and max bearing
  min_bear <- min(angles)
  max_bear <- max(angles)
  
  # if either are over 180, subtract 360
  if(min_bear > 180){
    min_bear <- min_bear - 360
  }
  if(max_bear > 180){
    max_bear <- max_bear - 360
  }
  
  
  # create lines to segment buffer
  max_end <- destPoint(p = col, b = max_bear, d = max_dist + 100000)
  min_end <- destPoint(p = col, b = min_bear, d = max_dist + 100000)
  max_line <- rbind(col, max_end) %>%
    as.lines(crs = "epsg:4326")
  min_line <- rbind(col, min_end) %>%
    as.lines(crs = "epsg:4326")
  lines <- aggregate(rbind(min_line, max_line))
  
  # erase lines from buffer and split buffer into segments
  lines <- buffer(lines, 0.01)
  buff <- erase(buff, lines)
  buff <- disagg(buff)
  buff$id <- 1:2
  
  # select segment containing tracks
  whichbuff <- extract(buff, tracks)
  buffno <- whichbuff %>%
    group_by(id) %>%
    summarise(n = n()) %>%
    filter(n == max(n)) %>%
    pull(id)
  buff <- buff %>% filter(id == buffno)
  
  
  # 4. Erase Any Land from Buffer
  
  # crop coast to extent of buffer
  crop_coast <- crop(coast, ext(buff))
  
  # erase land 
  buff <- erase(buff, crop_coast)
  
  # plot final boundary with tracks to check background sampling for each colony
  p1 <- ggplot() + geom_spatvector(data = buff, fill = "goldenrod1") +
    geom_spatvector(data = tracks, size = 0.5) + 
    geom_spatvector(data = crop(coast, ext(tracks))) +
    theme_bw()
  
  # export plot for checking
  ggsave(filename = paste0("output/background/plots/", this.species, "/", this.site, " ", this.stage, " boundary.png"), 
         plot = p1, width = 10, height = 10, units = "in", dpi = 300)
  
  
  # 5. Sample Background Points Relative to Distances from Colony
  
  # within each 5% interval of total distance to colony, what proportion of tracks fall into each band
  distance_props <- tracks$distances %>% cut(., breaks = seq(0, max_dist, length.out = 21)) %>%
    table() %>%
    prop.table() %>%
    data.frame() %>%
    mutate(band = 1:20)
  
  # loop over each band
  for(i in 1:20){
    
    # create buffer corresponding to 5% distance banding
    buff_band <- buffer(colony, max_dist * 0.05 * i)
    
    # only keep region that intersects with segment
    buff_band <- intersect(buff_band, buff)
    
    # if not the first band, erase the previous band
    if(i > 1){
      buff_band_erased <- erase(buff_band, erasure)
    } else {
      buff_band_erased <- buff_band
    }
    
    # bank this band for erasing in the next iteration
    erasure <- buff_band
    
    # if more than 20,000 track locations, sample nrow(tracks), otherwise sample 20,000
    if(nrow(tracks) > 20000){
      sample_size <- nrow(tracks)
    } else {
      sample_size <- 20000
    }
    
    # get number of points to sample from distance_props table
    n <- distance_props %>%
      filter(band == i) %>%
      pull(Freq) * sample_size
    n <- round(n)
    
    # skip if n = 0 
    if(n == 0){
      next
    }
    
    # sample within this band
    back_band <- spatSample(buff_band_erased, n) 
    
    # join to all other background points
    if(!exists("back")){
      back <- back_band
    } else {
      back <- bind_spat_rows(back, back_band)
    }
  }
  
  # assign dates and individual IDs proportionately from tracks
  back <- back %>%
    mutate(date = sample(tracks$date, n(), replace = T),
           individual_id = sample(tracks$individual_id, n(), replace = T))
  
  
  # 6. Extract Covariate Values
  
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
  
  # convert to dataframe
  back <- back %>%
    as.data.frame(geom = "XY")
  
  # 7. Export
  
  # save
  saveRDS(back, paste0("output/background/", this.species, "/", this.site, " ", this.stage, " background.rds"))
  
  # print completion
  print(paste0(this.species, " ", this.site, " ", this.stage, " complete"))
  
}

