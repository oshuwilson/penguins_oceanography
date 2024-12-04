#extract environmental variables for RSF
#ADD FSLE

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
}


# 1. Setup

#define species site and stage
this.species <- "ADPE"
this.site <- "Ardley Island"
this.stage <- "chick-rearing"

#load in tracks and background data
tracks <- readRDS(paste0("data/tracks/", this.species, "/", this.site, "/", this.stage, ".RDS"))
back <- readRDS(paste0("data/tracks/", this.species, "/", this.site, "/", this.stage, "_background.RDS"))

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


# 2.2 Dynamic Variables (add chl or wind?)

#load in dynamic_extract functions
source("code/functions/extraction_functions.R")

#sst 
data <- dynamic_extract("sst", data)

#mld
data <- dynamic_extract("mld", data)

#sal
data <- dynamic_extract("sal", data)

#ssh
data <- dynamic_extract("ssh", data)

#sic
data <- dynamic_extract("sic", data)
data$sic[is.na(data$sic)] <- 0 #SIC values of 0 print as NA in GLORYS

#curr
data <- dynamic_extract("uo", data)
data <- dynamic_extract("vo", data)
data$curr <- sqrt((data$uo^2) + (data$vo^2))


# 2.3 Oceanographic Variables

#front_freq
data <- dynamic_extract("front_freq", data)

#eddies
if(year(max(data$date)) > 2021){
  sub_data <- data %>% 
    filter(year(date) > 2021) %>%
    mutate(eddies = NA)
  data <- data %>% 
    filter(year(date) <= 2021)
  data <- dynamic_extract("eddies", data)
  data <- bind_spat_rows(data, sub_data)
} else{
  data <- dynamic_extract("eddies", data)
}

#dist2ice
data <- dynamic_extract("dist2ice", data)

#leads
leads <- rast("E:/Satellite_Data/static/leads/leads_resampled.nc")
data$leads <- extract(leads, data, ID=F)
rm(leads)


# 3. Format for export

#convert to data frame
data <- as.data.frame(data, geom = "XY")

#export
saveRDS(data, file = paste0("output/extractions/", this.species, "/", this.site, "_", this.stage, "_extracted.RDS"))

