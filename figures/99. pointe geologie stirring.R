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
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/ADPE/Pointe Geologie incubation tracks checked.rds")

# convert to terra
trax <- tracks %>%
  arrange(individual_id, date) %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# get dates for these tracks
dates <- sort(unique(as_date(trax$date)))

# read in salinity
sal <- rast("E:/Satellite_Data/daily/sal/sal_2017.nc")

# filter sal to example date
sal <- sal[[time(sal) == as_date("2017-12-13")]]

# crop to extent of tracks
e <- ext(trax) + c(0.5, 0.5, 0.5, 0.5)
sal <- crop(sal, e)

# truncate
sal <- clamp(sal, 31, 35)

# plot
ggplot() + geom_spatraster(data = sal) +
  facet_wrap(~lyr) +
  scale_fill_gradient()

# pick out examples
ex <- sal

# get depth for bathymetry
bathy <- getNOAA.bathy(lon1 = e[1], lon2 = e[2], lat1 = e[3], lat2 = e[4], resolution = 1)
bathy <- as.xyz(bathy) %>%
  rename(lon = V1, lat = V2, depth = V3) %>%
  rast()
plot(bathy)

# get contours
contours <- as.contour(bathy, levels = seq(-5000, 0, 1000))
crs(contours) <- "epsg:4326"
plot(contours, add = T)

# truncate values over 6 to 6
ex <- clamp(ex, 33, 35)

# get tracks for this year
year <- unique(year(time(ex)))
year_trax <- trax %>%
  filter(year(date) == year)

# create outline of land
land <- contours %>%
  filter(level == 0)

# plot individually
p1 <- ggplot() + geom_spatraster(data = ex) +
  geom_spatvector(data = contours, col = "black") +
  geom_spatvector(data = year_trax, col = "white", size = 1, alpha = 0.6) +
  geom_spatvector(data = land, col = "black", lwd = 1) +
  scale_fill_gradient2(midpoint = 33.8, mid = "grey80") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(134.2, 147.4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-67.3, -62.7))
p1
