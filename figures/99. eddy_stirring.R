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
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(rnaturalearthhires)
}


# 1. Chick-Rearing Macaronis from Fairy Point

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/MAPE/Fairy Point, Bird Island early chick-rearing tracks checked.rds")

# convert to terra
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# remove extra long track
trax <- trax %>%
  filter(individual_id != "MAPE-dtsetBirdLife751-H35-RAATD")

# limit tracks to high eddy year
high_trax <- trax %>%
  filter(year(date) == 2004)

# get dates for these tracks
high_dates <- sort(unique(as_date(high_trax$date)))

# limit tracks to within a week of eddy raster
rel_inds <- high_trax %>%
  filter(date > as_date("2004-01-23")) %>%
  pull(individual_id) %>%
  unique()
high_trax <- high_trax %>%
  filter(individual_id %in% rel_inds)

# keep individuals that forage at eddy periphery
high_trax <- high_trax %>%
  filter(!individual_id %in% c("MAPE-dtsetBirdLife751-H40-RAATD", "MAPE-dtsetBirdLife751-H42-RAATD"))

# read in nekton for high-eddy year
high_nek <- rast("E:/Satellite_Data/daily/epipelagic_nekton/epipelagic_nekton_2004.nc")

# filter nekton to these dates
high_nek <- high_nek[[time(high_nek) %in% high_dates]]

# crop to extent of tracks
e <- ext(-41, -37, -54.2, -52)
high_nek <- crop(high_nek, e)

# pick out examples
high_ex <- high_nek[[25]]

# get depth for bathymetry
bathy <- getNOAA.bathy(lon1 = e[1], lon2 = e[2], lat1 = e[3], lat2 = e[4], resolution = 1)
bathy <- as.xyz(bathy) %>%
  rename(lon = V1, lat = V2, depth = V3) %>%
  rast()

# truncate values over 6 to 6
high_ex <- clamp(high_ex, 0, 6)

# convert tracks to lines
high_lines <- as.lines(high_trax)

# create outline of South Georgia
sg <- ne_countries(scale = 10, returnclass = "sf") %>%
  vect()
sg <- crop(sg, e)


# make lines for each individual
inds <- unique(high_trax$individual_id)
for(ind in inds){
  
  # extract tracks
  ind_trax <- high_trax %>%
    filter(individual_id == ind)
  
  # create lines
  ind_lines <- as.lines(ind_trax)
  
  # join to all lines
  if(ind == inds[1]){
    all_lines <- ind_lines
  } else {
    all_lines <- c(all_lines, ind_lines)
  }
}
all_lines <- vect(all_lines)

# Fairy Point colony location
meta <- readRDS("data/metadata.RDS")
fp <- meta %>% 
  filter(deployment_site == "Fairy Point, Bird Island") %>%
  slice(1) %>%
  select(deployment_decimal_longitude, deployment_decimal_latitude)
fp <- vect(fp, geom = c("deployment_decimal_longitude", "deployment_decimal_latitude"), crs = "epsg:4326")

# plot individually
p1 <- ggplot() + geom_spatraster(data = high_ex) +
  geom_spatvector(data = all_lines, col = "grey30", lwd = 0.8) +
  geom_spatvector(data = sg, fill = "white") + 
  geom_spatvector(data = fp, size = 5, col = "black", fill = "red3", shape = 21) +
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "grey95", na.value = "grey95", midpoint = 1.5, 
                         name = bquote(atop(Nekton ~ Concentration, (g ~ m^-2)))) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(-40.295, -37.04)) +
  scale_y_continuous(expand = c(0, 0)) 
p1


# 2. Incubating Adelies from Pointe Geologie

# cleanup
rm(list=setdiff(ls(), "p1"))

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/ADPE/Pointe Geologie incubation tracks checked.rds")

# convert to terra
trax <- tracks %>%
  arrange(individual_id, date) %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# read in salinity
sal <- rast("E:/Satellite_Data/daily/sal/sal_2017.nc")

# filter sal to example date
sal <- sal[[time(sal) == as_date("2017-12-13")]]

# crop to extent of tracks
e <- ext(trax) + c(0.5, 0.5, 0.5, 0.5)
sal <- crop(sal, e)

# truncate
sal <- clamp(sal, 32.5, 35)

# get depth for bathymetry
bathy <- getNOAA.bathy(lon1 = e[1], lon2 = e[2], lat1 = e[3], lat2 = e[4], resolution = 1)
bathy <- as.xyz(bathy) %>%
  rename(lon = V1, lat = V2, depth = V3) %>%
  rast()

# get contours
contours <- as.contour(bathy, levels = seq(-5000, 0, 1000))
crs(contours) <- "epsg:4326"

# get tracks for this year
year <- unique(year(time(sal)))
year_trax <- trax %>%
  filter(year(date) == year)

# create outline of land
land <- contours %>%
  filter(level == 0)

# make lines for each individual
inds <- unique(year_trax$individual_id)
for(ind in inds){
  
  # extract tracks
  ind_trax <- year_trax %>%
    filter(individual_id == ind)
  
  # create lines
  ind_lines <- as.lines(ind_trax)
  
  # join to all lines
  if(ind == inds[1]){
    all_lines <- ind_lines
  } else {
    all_lines <- c(all_lines, ind_lines)
  }
}
all_lines <- vect(all_lines)

# plot individually
p2 <- ggplot() + geom_spatraster(data = sal) +
  #geom_spatvector(data = contours, col = "black") +
  geom_spatvector(data = all_lines, lwd = 0.8, col = "grey30") +
  #geom_spatvector(data = land, col = "black", lwd = 1) +
  scale_fill_gradient2(midpoint = 33.8, mid = "grey95", high = "darkred", low = "darkblue", name = "Salinity (PSU)",
                       na.value = "grey95") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(134.2, 147.4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-67.3, -62.7))
p2


# export
p1 <- p1 + theme(legend.position = c(1.12, 0.8))
ggsave("output/imagery/fig components/advection FP.png", p1, width = 11.2, height = 8, dpi = 300)

p2
ggsave("output/imagery/fig components/salinity PG.png", p2, width = 12, height = 8, dpi = 300)
