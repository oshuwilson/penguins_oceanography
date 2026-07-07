#-------------------------------------------------------------------------------
# Thinning ARS datasets
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
  library(GeoThinneR)
  library(ctmm)
  library(tidysdm)
}

# read in species site stage info to loop over
srs <- read.csv("data/tracks/species_site_stage_v2.csv")

# loop over each case study
for(i in 1:nrow(srs)){
  this.species <- srs$species[i]
  this.site <- srs$site[i]
  this.stage <- srs$stage[i]
  area <- srs$island[i]
  
  # print start 
  print(paste0("Starting ", this.species, " ", this.site, " ", this.stage))
  
  # load in hmm checked tracks
  tracks <- readRDS(paste0("output/extractions/", this.species, "/", this.site, " ", this.stage, " extracted.RDS")) %>%
    filter(pa == "presence")
  
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
  
  # set initial prop autocorr at 1
  prop.autocorr <- 1
  
  # set initial thinning interval as 0 hours
  thinterval <- 0
  
  # while prop.autocorr > 0.1
  while(prop.autocorr > 0.1){
    
    # add 1 hour to thinning interval
    thinterval <- thinterval + 1
    
    # limit to one point per every thinterval number of hours
    # round dates to nearest thinterval number of hours
    ars_thinned <- ars %>%
      mutate(date_round = round_date(date, unit = paste0(thinterval, " hours"))) %>%
      group_by(individual_id, date_round) %>%
      slice(1) %>%
      ungroup()
    
    # calculate autocorrelation for each individual
    inds <- unique(ars_thinned$individual_id)
    
    # start point of one for merging
    uno <- 1
    
    # loop
    for(ind in inds){
      
      # filter data
      ind.data <- ars_thinned %>% filter(individual_id == ind) %>%
        arrange(date)
      
      # make variables
      ind.data <- ind.data %>%
        mutate(lag.x = lag(x),
               lag.y = lag(y)) %>%
        mutate(dist = geosphere::distHaversine(cbind(x, y), cbind(lag.x, lag.y)))
      
      # if length of data is under 5, skip
      if(nrow(ind.data) < 5){
        uno <- uno + 1
        next
      }
      
      # calculate step length autocorrelation
      m1 <- acf(ind.data$dist[-1])
      
      # length of data
      n <- nrow(ind.data)
      
      # confidence intervals
      ci <- 1.96/sqrt(n)
      
      # extract autocorrelation for this individual
      autocorr <- m1$acf[2]
      
      # if autocorrelation is below confidence interval, okay
      if(abs(autocorr) < ci){
        autocorrelated <- 1
      }
      
      # create df
      df <- data.frame(id = ind,
                       autocorrelated = autocorrelated,
                       thinning_interval = thinterval)
      
      # bind to others
      if(ind == inds[uno]){
        all.df <- df
      } else {
        all.df <- rbind(all.df, df)
      }
    }
    
    # calculate proportion of individuals with autocorrelation
    prop.autocorr <- sum(all.df$autocorrelated)/nrow(all.df)
    
  }
  
  # use final thinning interval to thin all ARS data
  ars <- ars %>%
    mutate(date_round = round_date(date, unit = paste0(thinterval, " hours"))) %>%
    group_by(individual_id, date_round) %>%
    slice(1) %>%
    ungroup()
  
  # join back to background data
  back <- readRDS(paste0("output/extractions/", this.species, "/", this.site, " ", this.stage, " extracted.RDS")) %>%
    filter(pa == "absence")
  
  # final data for GAMMs
  data <- bind_rows(back, ars)
  
  # export dataframe
  saveRDS(data, 
          paste0("output/extractions/", this.species, "/", this.site, " ", this.stage, " thinned.RDS"))
  
}

