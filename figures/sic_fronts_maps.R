#-------------------------------------------
# Plot Eddy Specialisation by Sea Ice/Fronts
#-------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)

# read in species, site, stage info
srs <- read.csv("data/tracks/species_site_stage_v3.csv")

# remove unprocessed stages
srs <- srs %>% filter(eddy_importance != "")

# clean up south shetland and south orkney
srs <- srs %>%
  mutate(island = ifelse(island == "King George", "South Shetland", island)) %>%
  mutate(island = ifelse(island == "Signy", "South Orkney", island))

# summarise total and fraction of stages in each region that are specialised
fractions <- srs %>% 
  group_by(island) %>%
  summarise(total = n(),
            eddy_spec = sum(eddy_importance == "specialised"),
            fraction = sum(eddy_importance == "specialised")/n())

# read in metadata
meta <- readRDS("data/metadata.RDS")

# get colony location for first site from each stage
colony_locs <- srs %>% 
  group_by(island) %>%
  summarise(deployment_site = first(site)) 
colony_locs$deployment_site <- c(
  "Amanda Bay, Prydz Bay", "Cierva Cove, Antarctic Peninsula",
    "Atka Bay, Weddell Sea", "Auster Rookery, Mawson Coast",
    "Bechervaise Island, Mawson Coast", "Cape Bird, Ross Sea",
    "Cape Colbeck, Ross Sea", "Cape Crozier, Ross Sea", 
    "Cape Hallett, Ross Sea", "Baie du Marin, Crozet Islands",
    "Edmonson Point, Ross Sea", "Esperanza, Antarctic Peninsula",
    "Volunteer Beach, Falkland Islands", "Spit Bay, Heard Island",
    "Ratmanoff, Kerguelen Islands",
    "Sandy Bay, Macquarie Island", "Lutzow-Holm Bay, Mamejima Island",
    "Archway Bay, Marion Island", "Pointe Geologie, Terre Adelie",
    "Rothschild Island, Lazarev Bay", "Scullin Monolith, Mawson Coast",
    "Hound Bay, South Georgia",
    "Laurie Island, South Orkney", "Admiralty Bay, South Shetland Islands",
    "Taylor Glacier, Prince Charles Mountains", "Bahia Inutil, Tierra Del Fuego", 
    "Shirley Island, Windmill Islands")

# match with lat/lon from metadata
deps <- meta %>% 
  group_by(deployment_site) %>%
  summarise(lat = first(deployment_decimal_latitude), lon = first(deployment_decimal_longitude))
colony_locs <- colony_locs %>%
  left_join(deps)

# join locations with fractions
fractions <- fractions %>%
  left_join(colony_locs)

# convert to terra
frax <- fractions %>%
  vect(geom = c("lon", "lat"), crs = "epsg:4326") %>%
  project("epsg:3031")

# read in depth
depth <- rast("C:/Users/jcw2g17/OneDrive - University of Southampton/Documents/Predictor Data/depth_stereographic.RDS")
depth <- crop(depth, ext(-5790000, 5790000, -5790000, 5790000))

#classify depth into coarser bands
m <- c(0, 1000, 1,
       1000, 2000, 2,
       2000, 3000, 3,
       3000, 4000, 4,
       4000, 5000, 5,
       5000, 6000, 6,
       6000, 7000, 7,
       7000, 8000, 8,
       8000, 9000, 9)
m <- matrix(m, ncol = 3, byrow = TRUE)
depth <- classify(depth, m)

# read in fronts
fronts <- read_sf("exploration/colony summaries/south georgia/acc_fronts/shapefile/antarctic_circumpolar_current_fronts.shp")
fronts <- vect(fronts)

# remove subtropical front
fronts <- fronts %>% filter(NAME != "Subtropical Front (STF)")
fronts <- project(fronts, "epsg:3031")

# plot 
p1 <- ggplot() + geom_spatraster(data = depth) +
  geom_spatvector(data = fronts, lwd = 0.5) +
  geom_spatvector(data = frax, aes(col = fraction), size = 2.5) +
  theme_void() + 
  scale_fill_gradient(na.value = "white",
                      guide = "none",
                      high = "dodgerblue4",
                      low = "lightskyblue1") +
  scale_color_viridis_c(name = "Proportion")
p1

# export
ggsave("output/imagery/maps/front_map.png", p1, width = 8, height = 8, dpi = 300)

# read in monthly sic for 2000-2020 and limit to december-feb
sic <- rast(paste0("E:/Satellite_Data/monthly/sic/sic.nc"))
sic <- sic[[month(time(sic)) %in% c(1, 2, 12)]]
sic <- mean(sic, na.rm=T)
sic <- crop(sic, ext(-180, 180, -80, -50))
sic <- project(sic, "epsg:3031")
coast <- crop(coast, ext(sic))

p2 <- ggplot() + geom_spatraster(data = sic) +
  geom_spatvector(data = coast) +
  geom_spatvector(data = frax, aes(col = fraction), size = 2.5) +
  scale_fill_gradient(high = "white", low = "steelblue4", na.value = "steelblue4",
                      name = "Sea Ice Fraction") +
  scale_color_viridis_c(name = "Proportion") +
  theme_void()

# export
ggsave("output/imagery/maps/sic_map.png", p2, width = 8, height = 8, dpi = 300)
