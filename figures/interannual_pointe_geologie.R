#-------------------------
# Interannual Variability
#-------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(marmap)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(rnaturalearthhires)
}

# 2, Chick-Rearing Adelies from Pointe Geologie

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/ADPE/Pointe Geologie chick-rearing tracks checked.rds")

# convert to terra
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# round year to nearest year
trax <- trax %>%
  mutate(year = round_date(date, "year"))

# eddy foraging events
trax <- trax %>%
  mutate(edd_ars = ifelse(eddies >= 1 & state == "ARS" | eddies <= -1 & state == "ARS", 1, 0))
ggplot() + geom_spatvector(data = trax, aes(col = edd_ars))
ggplot(trax, aes(x = as.factor(year), y = edd_ars)) + geom_point()

# plot over mean SIC for each year with eddy foraging events
e <- ext(137, 142, -68, -64)
years <- c(2011, 2013, 2021, 2022)

for(this.year in years){
  last.year <- this.year-1
  # read in sic
  sic <- rast(paste0("E:/Satellite_Data/daily/sic/sic_", this.year, ".nc"))
  sic2 <- rast(paste0("E:/Satellite_Data/daily/sic/sic_", last.year, ".nc"))
  sic <- c(sic, sic2)
  
  # get dates for this year's tracks
  year_tracks <- trax %>%
    filter(year(year) == this.year)
  year_dates <- unique(as_date(year_tracks$date))
  
  # filter sic to these dates
  sic <- sic[[time(sic) %in% year_dates]]
  
  # crop sic
  sic <- crop(sic, e)
  
  # average sic
  sic <- mean(sic, na.rm = T)
  
  # plot
  p1 <- ggplot() + geom_spatraster(data = sic) +
    geom_spatvector(data = year_tracks, aes(col = edd_ars)) +
    scale_color_viridis_c() +
    ggtitle(this.year)
  print(p1)
}

# try for eddies instead
for(this.year in years[1:3]){
  last.year <- this.year - 1
  
  # read in eddies
  eddies <- rast(paste0("E:/Satellite_Data/daily/eddies/eddies_", this.year, ".nc"))
  eddies2 <- rast(paste0("E:/Satellite_Data/daily/eddies/eddies_", last.year, ".nc"))
  eddies <- c(eddies, eddies2)
  
  # get median date for this year's eddy ARS tracks
  year_tracks <- trax %>%
    filter(year(year) == this.year)
  edd_ars_year_tracks <- year_tracks %>%
    filter(edd_ars == 1)
  median_date <- median(as_date(edd_ars_year_tracks$date))
  
  # filter eddies to this date
  eddies <- eddies[[time(eddies) == median_date]]
  
  # crop eddies
  eddies <- crop(eddies, e)
  
  # plot
  p1 <- ggplot() + geom_spatraster(data = eddies) +
    geom_spatvector(data = year_tracks, aes(col = edd_ars)) +
    scale_color_viridis_c() +
    ggtitle(this.year)
  print(p1)
}
