#-------------------------------------------------------------------------------
# Match eddy amplitude, age, maturity, and intensity to locations
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(terra)
  library(gamm4)
  library(gratia)
  library(suncalc)
  library(tidysdm)
}

# read in species site stage info to loop over
srs <- read.csv("data/tracks/species_site_stage_v2.csv")

# choose in-depth case studies
for(i in c(8, 53, 54, 41)){
  
  this.species <- srs$species[i]
  this.site <- srs$site[i]
  this.stage <- srs$stage[i]
  area <- srs$island[i]
  
  # is rm_pack_ice yes 
  if(srs$rm_pack_ice[i] == "yes"){
    rm_pack_ice <- T
  } else {
    rm_pack_ice <- F
  }
  
  # load in hmm checked tracks
  tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))
  
  # if number of distinct individuals is less than 3, skip
  if(n_distinct(tracks$individual_id) < 3){
    next
  }
  
  # 1. Process Data
  
  # if fledglings, rename stage to post-breeding/pre-moult depending on species
  if(this.stage == "fledglings"){
    if(this.species == "ADPE"){
      og.stage <- "pre-moult"
    } else {
      og.stage <- "post-breeding"
    }
  } else {
    og.stage <- this.stage
  }
  
  #read in original tracks to get lat/lons and error info
  original <- readRDS(paste0("output/tracks/", this.species, "/", area, " ", og.stage, " tracks.RDS"))
  
  #append latitudes, longitudes, and errors to state tracks
  tracks <- tracks %>% 
    left_join(select(original, individual_id, date, lon, lat, 
                     longitude_se, latitude_se, 
                     lon_se_km, lat_se_km))
  
  #remove tracks with large error
  tracks <- tracks %>%
    filter((latitude_se < 0.05 & longitude_se < 0.125 | #greater allowance for longitude as this can be compressed at poles
              (lon_se_km < 5000 & lat_se_km < 5000)))
  
  #create column in date format for suncalc
  tracks <- tracks %>%
    rename(datetime = date) %>%
    mutate(date = as_date(datetime))
  
  #get dawn times
  tracks$dawn <- getSunlightTimes(data = tracks,
                                  keep = c("dawn"), tz = "UTC") %>%
    pull(dawn)
  
  #get dusk times
  tracks$dusk <- getSunlightTimes(data = tracks,
                                  keep = c("dusk"), tz = "UTC") %>%
    pull(dusk)
  
  #only keep points between sunrise and sunset
  tracks <- tracks %>%
    filter(datetime >= dawn & datetime <= dusk |
             is.na(dawn) & is.na(dusk))
  
  # filter to ARS only
  ars <- tracks %>% filter(state == "ARS")
  
  # resample non-eddies to 0
  ars <- ars %>%
    mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))
  
  # read in background samples
  back <- readRDS(paste0("output/background/", this.species, "/", this.site, " ", this.stage, " background.rds"))
  
  # resample non-eddies to 0
  back <- back %>%
    mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))
  
  # create binary presence/absence cols
  ars$pa <- 1
  back$pa <- 0
  
  # select key variables
  if(this.species %in% c("ADPE", "CHPE", "EMPE")){
    columns <- c("individual_id", "date", "ed2", "depth", "sic", "uo", "vo", "pa", "x", "y")
  } else {
    columns <- c("individual_id", "date", "ed2", "depth", "uo", "vo", "pa", "x", "y")
  }
  
  ars <- ars %>% 
    select(all_of(columns))
  back <- back %>%
    select(all_of(columns)) %>%
    mutate(date = as_date(date))
  
  # reproject ars to epsg:4326
  ars <- ars %>%
    vect(geom = c("x", "y"), crs = "epsg:6932") %>%
    project("epsg:4326") %>%
    as.data.frame(geom = "XY")
  
  # join datasets together
  data <- bind_rows(ars, back)
  
  # calculate EKE
  data <- data %>%
    mutate(eke = 0.5 * (uo^2 + vo^2))
  
  # remove sea ice concentrations above 10%
  if(rm_pack_ice == T){
    data <- data %>%
      filter(sic < 0.1)
  }
  
  
  # 2. Match to eddy IDs
  
  # read in eddy data
  cyc <- readRDS("data/eddies/cyclonic_attributes.RDS")
  anti <- readRDS("data/eddies/anticyclonic_attributes.RDS")
  
  # get unique dates in the data
  dates <- unique(as_date(data$date)) %>% sort()
  
  # limit eddy data to these dates
  cyc <- cyc %>%
    filter(time %in% dates)
  anti <- anti %>%
    filter(time %in% dates)
  
  # get lat/lon limits of data
  lon_min <- min(data$x) - 5
  lon_max <- max(data$x) + 5
  lat_min <- min(data$y) - 5
  lat_max <- max(data$y) + 5
  
  # normalise eddy longitudes
  eddies <- bind_rows(cyc %>% mutate(eddy_type = "cyclone"), 
                      anti %>% mutate(eddy_type = "anticyclone"))
  eddies <- eddies %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon))
  
  # limit eddy data to these lat/lon limits
  eddies <- eddies %>%
    filter(lat >= lat_min & lat <= lat_max &
             lon >= lon_min & lon <= lon_max)
  
  # normalise eddy age to get maturity
  eddies <- eddies %>%
    group_by(id) %>%
    mutate(eddy_maturity = observation_number / max(observation_number)) %>%
    ungroup()
  
  # create empty dataframe to store matched data
  data_eddies <- data.frame()
  
  # match data to eddies by date and location
  for(j in 1:length(dates)){
    this.date <- dates[j]
    
    # limit all eddies to this date
    eddies_sub <- eddies %>%
      filter(time == this.date)
    
    # limit locations to this date
    data_sub <- data %>%
      filter(date == this.date)
    
    # if nrow(data_sub) == 0, skip to next date
    if(nrow(data_sub) == 0){
      next
    }
    
    # create empty data frame
    eddy_ats <- data.frame()
    
    # for each location
    for(k in 1:nrow(data_sub)){
      loc <- data_sub[k,]
      
      # get lon/lat of location
      loc_lon <- loc$x
      loc_lat <- loc$y
      
      # limit eddies to those within 5 degrees of location
      eddies_loc <- eddies_sub %>%
        filter(lat >= loc_lat - 5 & lat <= loc_lat + 5 &
                 lon >= loc_lon - 5 & lon <= loc_lon + 5)
      
      # if nrow(eddies_loc) == 0, skip to next location
      if(nrow(eddies_loc) == 0){
        next
      }
      
      # calculate distance to each eddy
      eddies_loc$distances <- geodist::geodist(eddies_loc %>% select(lon, lat), 
                                               loc %>% select(x, y), 
                                               measure = "haversine")[,1]
      
      # normalise distances by eddy radius
      eddies_loc <- eddies_loc %>%
        mutate(norm_distance = distances / radius)
      
      # get minimum distance 
      eddies_loc <- eddies_loc %>% filter(norm_distance == min(norm_distance))
      
      # if minimum distance >2, skip to next location
      if(eddies_loc$norm_distance > 2){
        next
      } else{
        loc <- loc %>%
          mutate(eddy_id = eddies_loc$id,
                 eddy_type = eddies_loc$eddy_type,
                 eddy_amplitude = eddies_loc$amplitude,
                 eddy_radius = eddies_loc$radius,
                 eddy_age = eddies_loc$observation_number,
                 eddy_maturity = eddies_loc$eddy_maturity,
                 eddy_distance = eddies_loc$norm_distance)
        eddy_ats <- bind_rows(eddy_ats, loc)
      }
      
    }
    
    # if no eddy dataframe for this date, skip to next date
    if(nrow(eddy_ats) == 0){
      next
    } else {
      data_eddies <- bind_rows(data_eddies, eddy_ats)
    }
    
    # print progress
    print(paste0(j, "/", length(dates)))
    
  }
  
  print(ggplot(data_eddies) +
          geom_density(aes(x = eddy_amplitude, fill = as.factor(pa)), alpha = .5) +
          theme_bw())
  print(ggplot(data_eddies) +
          geom_density(aes(x = eddy_radius, fill = as.factor(pa)), alpha = .5) +
          theme_bw())
  print(ggplot(data_eddies) +
          geom_density(aes(x = eddy_age, fill = as.factor(pa)), alpha = .5) +
          theme_bw())  
  print(ggplot(data_eddies) +
          geom_density(aes(x = eddy_maturity, fill = as.factor(pa)), alpha = .5) +
          theme_bw())
  
  # export attributes
  saveRDS(data_eddies, paste0("output/eddy attributes/", this.species, " ", this.site, " ", this.stage, " eddy attributes.rds"))
  
}
