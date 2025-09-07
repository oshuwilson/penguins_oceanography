#-------------------------------------
# Plot South Georgia Front System
#-------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(sf)
  library(marmap)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(rnaturalearthhires)
}

# read in fronts
fronts <- read_sf("exploration/colony summaries/south georgia/sg_polar_fronts.shp")
fronts <- vect(fronts)
fronts <- project(fronts, "epsg:4326")

# load coastline
coast <- ne_countries(scale = 10, returnclass = "sf") %>%
  vect()

# create bounding box for south georgia
e <- ext(-43, -30, -57, -51)

# crop layers to extent
coast <- crop(coast, e)
fronts <- crop(fronts, e)

# get depth data
bathy <- getNOAA.bathy(lon1 = e[1], lon2 = e[2], lat1 = e[3], lat2 = e[4], resolution = 1)
bathy <- as.xyz(bathy) %>%
  rename(lon = V1, lat = V2, depth = V3) %>%
  rast()

# remove Polar Front - hardly in frame
fronts <- fronts %>%
  filter(!name %in% c("Polar Front (PF)", "Southern boundary of the Antarctic Circumpolar Current"))

# rename other fronts
fronts <- fronts %>%
  mutate(name = case_when(
    name == "Southern Antarctic Circumpolar Current Front (sACCf)" ~ "sACCf"
  ))

# read in Macaroni tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/MAPE/Fairy Point, Bird Island early chick-rearing tracks checked.RDS")

# convert to terra
trax <- vect(tracks, geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# remove one anomalous track
trax <- trax %>%
  filter(individual_id != "MAPE-dtsetBirdLife751-H35-RAATD")

# convert trax to lines one individual at a time
inds <- unique(trax$individual_id)
for(ind in inds){
  
  # get individual track
  trax_ind <- trax %>%
    filter(individual_id == ind)
  
  # convert to lines
  trax_ind <- trax_ind %>%
    as.lines()
  
  # join to other lines
  if(ind == inds[1]){
    trax_lines <- trax_ind
  } else {
    trax_lines <- c(trax_lines, trax_ind)
  }
}
trax_lines <- vect(trax_lines)

# Fairy Point colony location
meta <- readRDS("data/metadata.RDS")
fp <- meta %>% 
  filter(deployment_site == "Fairy Point, Bird Island") %>%
  slice(1) %>%
  select(deployment_decimal_longitude, deployment_decimal_latitude)
fp <- vect(fp, geom = c("deployment_decimal_longitude", "deployment_decimal_latitude"), crs = "epsg:4326")

# inset bounding box
bbox <- ext(-41, -37, -54.2, -52)
bbox <- as.polygons(bbox)
crs(bbox) <- "epsg:4326"

# plot
p1 <- ggplot() + 
  geom_spatraster(data = bathy) +
  geom_spatvector(data = coast, fill = "white") + 
  geom_spatvector(data = fronts, aes(linetype = name), lwd = 1.5) +
  scale_fill_gradient(limits = c(-8000, 0), name = "Depth (m)", 
                      high = "#eaeff3", low = "steelblue4", na.value = "#eaeff3") +
  geom_spatvector(data = fp, size = 3, col = "black", fill = "red3", shape = 21) +
  geom_spatvector(data = bbox, fill = NA, col = "grey40", lwd = 1) +
  scale_linetype(name = "", labels = "sACCf Position") +
  theme_classic() +
  theme(legend.spacing = unit(1, "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) 
p1

# export
ggsave("output/imagery/fig components/fairy_point_background.png", p1, width = 11, height = 8, units = "in", dpi = 300)
