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

# read in temperature at depths
temp <- rast("exploration/colony summaries/crozet/thetao2.nc")

# crop to extent
e <- ext(40, 60, -55, -44)
temp2 <- crop(temp, e)

# clamp to prevent high values in northeast corner warping it
temp2 <- clamp(temp2, 0, 7)

# load in land
land <- ne_countries(scale = 10, returnclass = "sv")
land <- crop(land, e)

# create bounding box
bbox <- ext(50.205, 52.795, -48.378, -46.125)
bbox <- as.polygons(bbox)
crs(bbox) <- "epsg:4326"

# get colony location
meta <- readRDS("data/metadata.RDS")
pos <- meta %>% 
  filter(deployment_site == "Ile de la Possession, Crozet") %>%
  select(deployment_decimal_longitude, deployment_decimal_latitude) %>%
  vect(geom = c("deployment_decimal_longitude", "deployment_decimal_latitude"), crs = "epsg:4326")


# plot all together
p1 <- ggplot() + 
  geom_spatraster(data = temp2) + 
  geom_spatvector(data = land, fill = "white") +
  geom_spatvector(data = pos, size = 6, col = "black", fill = "red3", shape = 21) +
  geom_spatvector(data = bbox, fill = NA, col = "grey70", lwd = 1.5) +
  scale_fill_viridis_c(name = "Temperature \nat 155m (°C)") +
  theme_classic() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) 
p1

ggsave("output/imagery/fig components/crozet background.png", p1,
       width = 12, height = 10, dpi = 300)
