#---------------------------
# ARS Prevalence GAMMs
#---------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(mgcv)
  library(gamm4)
  library(ggside)
  library(patchwork)
  library(gratia)
  library(foreach)
  library(suncalc)
}

#function to fit GAMM using covariates with appropriate unique values
pred_gam <- function(tracks){
  gamm4::gamm4(
    as.formula(
      paste0(
        "bin_state ~ s(", 
        vars %>% paste0(collapse = ", bs = 'ts', k=5) + s("),
        ", bs = 'ts', k=5)")
    ), 
    random = ~(1|individual_id), family=binomial, data=tracks)
}

# 1. Data Preparation

# define species
this.species <- "KIPE"
longname <- "King Penguin"

# read in species, region, and stage info
srs <- read.csv("data/tracks/species_site_stage_v2.csv")
srs <- srs %>%
  filter(species == this.species)

# loop over stages
stages <- unique(srs$stage)
for(this.stage in stages){
  
  # loop over regions
  regions <- srs %>%
    filter(stage == this.stage) %>%
    pull(site)
  
  for(this.site in regions){
    try({
    
    # find larger area to get original tracks
    area <- srs %>% 
      filter(site == this.site, stage == this.stage) %>%
      pull(island)
      
    # load the tracks with states
    tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))
    
    # if fewer than 3 individuals, skip to next site
    if(length(unique(tracks$individual_id)) < 3){
      next
    }
    
    # change state to binomial variable
    tracks <- tracks %>%
      mutate(bin_state = if_else(state == "ARS", 1, 0))
    
    #resample eddy 0s to values between -0.99 and 0.99 - to avoid misleading trends between -1 and 1
    tracks <- tracks %>%
      group_by(row_number()) %>%
      mutate(ed2 = ifelse(eddies == 0, runif(n(), -0.99, 0.99), eddies)) %>%
      ungroup()
    
    #read in original tracks to get lat/lons and error info
    original <- readRDS(paste0("output/tracks/", this.species, "/", area, " ", this.stage, " tracks.RDS"))
    
    #append latitudes, longitudes, and errors to state tracks
    tracks <- tracks %>% 
      left_join(select(original, individual_id, date, lon, lat, 
                       longitude_se, latitude_se, 
                       lon_se_km, lat_se_km))
    
    #remove tracks with large error
    tracks <- tracks %>%
      filter((latitude_se < 0.133 & longitude_se < 0.133) |
               (lon_se_km < 9000 & lat_se_km < 9000))
    
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
    
    #set seed for reproducibility
    set.seed(777)
    
    
    # 2. Fit GAMMs
    
    #list of all possible covariates
    allvars <- c("ed2", "curr", "depth")
    
    #if ADPE or EMPE, account for sea ice foraging
    if(this.species %in% c("ADPE", "EMPE")){
      allvars <- c("ed2", "curr", "depth", "sic")
    }
    
    #find covariates with too little variation for modelling (fewer than 6 unique values)
    lowvar <- tracks %>%
      select(all_of(allvars)) %>% 
      summarise_all(n_distinct) %>%
      pivot_longer(cols = everything(), names_to = "covariate", values_to = "unique") %>%
      filter(unique < 10) %>%
      pull(covariate)
    
    #remove these covariates
    vars <- allvars[!allvars %in% lowvar]
    
    #remove covariates with unsuitable variation between individuals (i.e. only one individual shows variation)
    #typically only sea ice concentration
    if(this.species %in% c("ADPE", "EMPE")){
      smallvar <- tracks %>%
        select(sic, individual_id) %>%
        group_by(individual_id) %>%
        summarise_all(n_distinct)
      
      ss <- smallvar %>% 
        filter(sic == 1) %>%
        nrow()
      if(ss >= nrow(smallvar) - 1){
        vars <- vars[vars != "sic"]
      }
    }
    
    #fit GAMM using custom function
    m1 <- pred_gam(tracks)
    
    # get p-values
    m1sum <- as.list(summary(m1$gam))
    pvals <- m1sum$s.pv
    
    # create table of covariate names and p-values
    cov_names <- vars
    pscores <- data.frame(cov_names, pvals)
    
    # get smooths
    sm <- smooth_estimates(m1$gam) %>%
      add_confint()
    
    # add constant
    constant <- coef(m1$gam)[1]
    sm <- sm %>%
      mutate(.estimate = .estimate + constant,
             .lower_ci = .lower_ci + constant,
             .upper_ci = .upper_ci + constant)
    
    
    # transform the smooths
    backtrans <- inv_link(m1$gam)
    sm <- sm %>%
      mutate(.estimate = backtrans(.estimate),
             .lower_ci = backtrans(.lower_ci),
             .upper_ci = backtrans(.upper_ci))
    
    # add site and stage to dataframe
    sm <- sm %>%
      mutate(site = this.site,
             stage = this.stage)
    
    # add NAs for any covariates not used in the GAMM
    missing <- allvars[!allvars %in% cov_names]
    for(missingvar in missing){
      sm <- sm %>%
        mutate(!!missingvar := NA)
    }
    
    #return the sm dataframe
    if(this.stage == stages[1] & this.site == regions[1]){
      smooths <- sm
    } else {
      smooths <- rbind(smooths, sm)
    }
    
    #return the p-values as a dataframe
    pscores <- pscores %>%
      mutate(site = this.site,
             stage = this.stage)
    if(this.stage == stages[1] & this.site == regions[1]){
      pvalues <- pscores
    } else {
      pvalues <- rbind(pvalues, pscores)
    }
    })
  }
}

#export smooths
saveRDS(smooths, paste0("output/gamms/smooths/", this.species, "_smooths.rds"))

#export 
saveRDS(pvalues, paste0("output/gamms/pvalues/", this.species, "_pvalues.rds"))
