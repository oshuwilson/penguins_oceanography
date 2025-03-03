#----------------------------------------
# Plots Showing Examples of Eddy Stirring
#----------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
}


# 1. Pointe Geologie Incubating Adelies

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/ADPE/Pointe Geologie incubation tracks checked.rds")

# try specific years - 2017 looks good
year <- 2017

# filter tracks to that year
tracks <- tracks %>% 
  filter(year(date) == year)

# read in salinity for that year
sal <- rast(paste0("E:/Satellite_Data/daily/sal/sal_", year, ".nc"))
sal <- rast(paste0("E:/Satellite_Data/daily/sic/sic_", year, ".nc"))


# limit salinity to dates in tracks
dates <- unique(as_date(tracks$date))
sal <- sal[[time(sal) %in% dates]]

# convert tracks to terra
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")
e <- ext(trax)
e

# find track IDs that use eddies
# sort(unique(trax$individual_id))
# trax <- trax %>%
#   filter(individual_id %in% sort(unique(trax$individual_id))[1:4])
ggplot() + geom_spatvector(data = trax, aes(col = individual_id))

# probably A12_F1_20171211_CT or A1_F1_20171204_CT
# start with A12_F1_20171211_CT
ind <- trax %>%
  filter(individual_id == "A12_F1_20171211_CT")

# find median date of this track
med_date <- median(as_date(ind$date)) - days(10)
med_date

# limit salinity to this date 
ind_sal <- sal[[time(sal) == med_date]]

# crop salinity to this track's extent
e <- ext(ind) + c(2, 2, 2, 2)
ind_sal <- crop(ind_sal, e)

# plot
ggplot() + geom_spatraster(data = ind_sal) +
  geom_spatvector(data = trax) 


sal <- crop(sal, e)
plot(sal)


# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/ADPE/Pointe Geologie incubation tracks checked.rds")

# try specific years - 2017 looks good
year <- 2011

# filter tracks to that year
tracks <- tracks %>% 
  filter(year(date) == year)

# read in salinity for that year
sal <- rast(paste0("E:/Satellite_Data/daily/sal/sal_", year, ".nc"))
sal <- rast(paste0("E:/Satellite_Data/daily/chl/resampled/chl_", year, "_resampled.nc"))


# limit salinity to dates in tracks
dates <- unique(as_date(tracks$date))
sal <- sal[[time(sal) %in% dates]]

# convert tracks to terra
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")
e <- ext(trax)
e

# find track IDs that use eddies
# sort(unique(trax$individual_id))
# trax <- trax %>%
#   filter(individual_id %in% sort(unique(trax$individual_id))[1:4])
ggplot() + geom_spatvector(data = trax, aes(col = individual_id))

# find median date of this track
med_date <- median(as_date(trax$date))
med_date

# limit salinity to this date 
ind_sal <- sal[[time(sal) == med_date]]

# crop salinity to this track's extent
e <- ext(trax) + c(2, 2, 2, 2)
ind_sal <- crop(ind_sal, e)

# plot
ggplot() + geom_spatraster(data = ind_sal) +
  geom_spatvector(data = trax) +
  scale_fill_viridis_c(trans = "log")


sal <- crop(sal, e)
ggplot() + geom_spatraster(data = sal) +
  scale_fill_viridis_c(trans = "log") +
  facet_wrap(~lyr)


# TRY MAPES