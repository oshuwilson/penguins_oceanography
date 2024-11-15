# extract covariates for observed and available steps for SSFs

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
}


# 1. Setup

#define species site and stage
this.species <- "KIPE"
this.site <- "Crozet"
this.stage <- "incubation"

#read in steps
steps <- readRDS(paste0("output/ssf/available steps/", this.species, "/", this.site, "_", this.stage, "_available_steps.RDS"))

#format date column
steps <- steps %>%
  mutate(date = t2_)

#create SpatVector for steps
steps <- vect(steps,
              geom=c("x2_", "y2_"),
              crs="epsg:6932")


# 2. Extract environmental variables

# 2.1 Static Variables

#depth
depth <- rast("E:/Satellite_Data/static/depth/depth.nc")

#project steps to same CRS as depth
steps <- project(steps, crs(depth))

#extract
steps$depth <- extract(depth, steps, ID=F)

#remove rows where depth is NA - will be NA for every GLORYS variable
steps <- steps %>% drop_na(depth)

#slope
slope <- rast("E:/Satellite_Data/static/slope/slope.nc")
steps$slope <- extract(slope, steps, ID=F)

#dshelf
dshelf <- rast("E:/Satellite_Data/static/dshelf/dshelf_resampled.nc")
steps$dshelf <- extract(dshelf, steps, ID=F)

#cleanup static
rm(depth, slope, dshelf)

# 2.2 Dynamic Variables (add chl or wind?)

#load in dynamic_extract functions
source("code/functions/extraction_functions.R")

#sst 
steps <- dynamic_extract("sst", steps)

#mld
steps <- dynamic_extract("mld", steps)

#sal
steps <- dynamic_extract("sal", steps)

#ssh
steps <- dynamic_extract("ssh", steps)

#sic
steps <- dynamic_extract("sic", steps)
steps$sic[is.na(steps$sic)] <- 0 #SIC values of 0 print as NA in GLORYS

#curr
steps <- dynamic_extract("uo", steps)
steps <- dynamic_extract("vo", steps)
steps$curr <- sqrt((steps$uo^2) + (steps$vo^2))

#front_freq
steps <- dynamic_extract("front_freq", steps)

#eddies
steps <- dynamic_extract("eddies", steps)

#dist2ice
steps <- dynamic_extract("dist2ice", steps)

#leads
leads <- rast("E:/Satellite_Data/static/leads/leads_resampled.nc")
steps$leads <- extract(leads, steps, ID=F)
rm(leads)


# 3. Export

#reproject steps to planar view
steps <- project(steps, "epsg:6932")

#extract data frame
data <- as.data.frame(steps, geom = "XY")

#format
data <- data %>% 
  rename(x2_ = x, y2_ = y) %>%
  select(-date, -uo, -vo)

#export
saveRDS(data, paste0("output/ssf/extracted/", this.species, "/", this.site, "_", this.stage, "_steps_extracted.RDS"))

        