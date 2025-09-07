#-------------------------------------
# Plot Pointe Geologie Sea Ice Eddies
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

# load coastline
coast <- ne_countries(scale = 10, returnclass = "sf") %>%
  vect()

# load sea ice for set day
sic <- rast("E:/Satellite_Data/daily/sic/sic_2017.nc")
sic <- sic[[time(sic) == as_date("2017-12-13")]]

# create bounding box for pointe geologie
e <- ext(125, 160, -70, -60)

# crop layers to extent
coast <- crop(coast, e)
sic <- crop(sic, e)
plot(sic)
plot(coast, add = T)

# read in Adelie tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/ADPE/Pointe Geologie incubation tracks checked.RDS")

# convert to terra
trax <- vect(tracks, geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# only keep tracks active at this time
active_inds <- trax %>% 
  filter(as_date(date) == as_date("2017-12-13")) %>%
  pull(individual_id) %>%
  unique()
trax <- trax %>%
  filter(individual_id %in% active_inds)
plot(trax, add = T, pch = ".")

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
plot(trax_lines, add = T, col = "red")

# Pointe Geologie colony location
meta <- readRDS("data/metadata.RDS")
pg <- meta %>% 
  filter(deployment_site == "Pointe Geologie, Terre Adelie") %>%
  slice(1) %>%
  select(deployment_decimal_longitude, deployment_decimal_latitude)
pg <- vect(pg, geom = c("deployment_decimal_longitude", "deployment_decimal_latitude"), crs = "epsg:4326")

# polygon for cropping ADD coastfile
bbox2 <- as.polygons(e)
crs(bbox2) <- "epsg:4326"
bbox2 <- project(bbox2, "epsg:3031")

# read in coast and project
coast <- sf::read_sf("~/OneDrive - University of Southampton/Documents/Predictor Data/coast/add_coastline_medium_res_polygon_v7_9.shp")
coast <- vect(coast)
coast <- crop(coast, bbox2)
coast <- project(coast, "epsg:4326")
plot(coast)

# new bbox 
bbox <- as.polygons(ext(136, 150, -68, -62))
crs(bbox) <- "epsg:4326"

# read in depth
depth <- rast("E:/Satellite_Data/static/depth/depth.nc")

# crop depth
depth <- crop(depth, e)

# make all sic values of NA into 0
values(sic)[is.na(values(sic))] <- 0
plot(sic)

# make mask of NA values
sic <- mask(sic, depth)

# plot
p1 <- ggplot() + 
  geom_spatraster(data = sic) +
  geom_spatvector(data = coast, fill = "white") + 
  #geom_spatvector(data = shelf) +
  scale_fill_gradient(limits = c(0, 1), name = "Sea Ice Fraction", 
                      high = "white", low = "steelblue4", na.value = "grey60") +
  geom_spatvector(data = bbox, fill = NA, col = "grey60", lwd = 1) +
  geom_spatvector(data = pg, size = 5, col = "black", fill = "red3", shape = 21) +
  theme_classic() +
  theme(legend.spacing = unit(1, "cm")) +
  scale_y_continuous(limits = c(-69.9, -60), expand = c(0,0)) +
  scale_x_continuous(limits = c(125.85, 158), expand = c(0,0)) 
p1

# crop everything to bbox
sic2 <- crop(sic, bbox)
coast2 <- crop(coast, bbox)

# plot tracks
p2 <- ggplot() +
  geom_spatraster(data = sic2) +
  geom_spatvector(data = coast2, fill = "white") + 
  geom_spatvector(data = trax_lines, col = "black", lwd = 0.8) +
  geom_spatvector(data = pg, size = 8, col = "black", fill = "red3", shape = 21) +
  scale_fill_gradient(limits = c(0, 1), name = "Sea Ice Fraction", 
                      high = "white", low = "steelblue4", na.value = "grey60") +
  theme_classic() +
  scale_x_continuous(expand = c(0,0), limits = c(136.2, 149.8)) +
  scale_y_continuous(expand = c(0,0), limits = c(-67.88, -62.2))
p2

# get shelf break
# shelf <- sf::read_sf("output/imagery/fig components/antarctic_shelf_break/Shelf_break_Antarctica_DAmblas.shp")
# shelf <- vect(shelf)
# bbox2 <- project(bbox2, crs(shelf))
# shelf <- crop(shelf, bbox2)
# shelf <- project(shelf, "epsg:4326")

# export
ggsave("output/imagery/fig components/pointe_geologie_background.png", p1, width = 11, height = 8, units = "in", dpi = 300)
ggsave("output/imagery/fig components/pointe_geologie_tracks.png", p2, width = 8, height = 8, units = "in", dpi = 300)
