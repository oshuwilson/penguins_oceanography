#----------------------------------------
# Plots Showing Examples of Eddy Stirring
#----------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(marmap)
  library(cowplot)
}

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/MAPE/Fairy Point, Bird Island early chick-rearing tracks checked.rds")

# convert to terra
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# remove extra long track
trax <- trax %>%
  filter(individual_id != "MAPE-dtsetBirdLife751-H35-RAATD")

# limit tracks to high and low eddy years
high_trax <- trax %>%
  filter(year(date) == 2004)
low_trax <- trax %>%
  filter(year(date) == 2012)

# get dates for these tracks
high_dates <- sort(unique(as_date(high_trax$date)))
low_dates <- sort(unique(as_date(low_trax$date)))

# read in nekton for high-eddy and low-eddy year
high_nek <- rast("E:/Satellite_Data/daily/epipelagic_nekton/epipelagic_nekton_2004.nc")
low_nek <- rast("E:/Satellite_Data/daily/epipelagic_nekton/epipelagic_nekton_2012.nc")

# filter nekton to these dates
high_nek <- high_nek[[time(high_nek) %in% high_dates]]
low_nek <- low_nek[[time(low_nek) %in% low_dates]]

# crop to extent of tracks
e <- ext(trax) + c(0.5, 0.5, 0.5, 0.5)
high_nek <- crop(high_nek, e)
low_nek <- crop(low_nek, e)

# plot
ggplot() + geom_spatraster(data = high_nek) +
  #geom_spatvector(data = high_trax, size = 0.1) +
  facet_wrap(~lyr) +
  scale_fill_viridis_c(limits = c(0, 12))
ggplot() + geom_spatraster(data = low_nek) +
  #geom_spatvector(data = low_trax, size = 0.1) +
  facet_wrap(~lyr) +
  scale_fill_viridis_c(limits = c(0, 12))

# pick out examples
high_ex <- high_nek[[25]]
low_ex <- low_nek[[16]]

# get depth for bathymetry
bathy <- getNOAA.bathy(lon1 = e[1], lon2 = e[2], lat1 = e[3], lat2 = e[4], resolution = 1)
bathy <- as.xyz(bathy) %>%
  rename(lon = V1, lat = V2, depth = V3) %>%
  rast()
plot(bathy)

# get contours
contours <- as.contour(bathy, levels = seq(-5000, 0, 1000))
crs(contours) <- crs(trax)

# truncate values over 6 to 6
high_ex <- clamp(high_ex, 0, 6)
low_ex <- clamp(low_ex, 0, 6)

# convert tracks to lines
high_lines <- as.lines(high_trax)
low_lines <- as.lines(low_trax)

# create outline of South Georgia
sg <- contours %>%
  filter(level == 0)

# plot individually
p1 <- ggplot() + geom_spatraster(data = high_ex) +
  geom_spatvector(data = contours, col = "black") +
  geom_spatvector(data = high_lines, col = "white", lwd = 1, alpha = 0.6) + 
  scale_fill_viridis_c(limits = c(0, 6), 
                       name = bquote(Nekton ~ Concentration ~ (g ~ m^-2))) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) 
p1

p2 <- ggplot() + geom_spatraster(data = low_ex) +
  geom_spatvector(data = contours, col = "black") +
  geom_spatvector(data = low_lines, col = "white", lwd = 1, alpha = 0.6) + 
  scale_fill_viridis_c(limits = c(0, 6), 
                       name = bquote(Nekton ~ Concentration ~ (g ~ m^-2))) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) 
p2

p1 <- ggplot() + geom_spatraster(data = high_ex) +
  geom_spatvector(data = contours, col = "black") +
  geom_spatvector(data = high_lines, col = "white", lwd = 2, alpha = 0.6) +
  geom_spatvector(data = sg, col = "black", lwd = 1) +
  scale_fill_gradient(limits = c(0, 6), low = "grey90", na.value = "grey90", high = "pink4",
                       name = bquote(Nekton ~ Concentration ~ (g ~ m^-2))) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(-40.295, -37.04)) +
  scale_y_continuous(expand = c(0, 0)) 
p1

