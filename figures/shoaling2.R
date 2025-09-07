rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(marmap)
}

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/KIPE/Ile de la Possession, Crozet chick-rearing tracks checked.rds")
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# filter to year with most tracks
sub_trax <- trax %>%
  filter(year(date) == 2013)

# create lines for each individual
inds <- unique(sub_trax$individual_id)
for(ind in inds){
  
  ind_trax <- sub_trax %>%
    filter(individual_id == ind)
  ind_trax <- as.lines(ind_trax)
  plot(ind_trax)
  if(ind == inds[1]){
    trax_lines <- ind_trax
  } else {
    trax_lines <- c(trax_lines, ind_trax)
  }
}
trax_lines <- vect(trax_lines)

# read in temperature at depths
temp <- rast("exploration/colony summaries/crozet/thetao.nc")

# limit temporally
temp <- temp[[time(temp) == as_date("2013-02-05")]]

# limit to 155m depth
temp <- temp[[25]]

# get colony location
meta <- readRDS("data/metadata.RDS")
unique(meta$deployment_site)
pos <- data.frame(x = 51.87, y = -46.43) %>%
  vect(geom = c("x", "y"), crs = "epsg:4326")
pos  

# load in land
land <- ne_countries(scale = 10, returnclass = "sv")
land <- crop(land, ext(50, 53, -49, -46))
plot(land)

# plot
p1 <- ggplot() + geom_spatraster(data = temp) +
  geom_spatvector(data = land, fill = "white") +
  geom_spatvector(data = trax_lines, lwd=1) +
  scale_fill_viridis_c(name = "Temperature \nat 155m (°C)") +
  theme_classic() +
  scale_y_continuous(limits = c(-48.378, -46.125), expand = c(0,0)) +
  scale_x_continuous(limits = c(50.205, 52.795), expand = c(0,0)) +
  geom_spatvector(data = pos, size = 6, col = "black", fill = "red3", shape = 21) 
p1

ggsave("output/imagery/fig components/thermocline.png", p1,
       width = 9, height = 11, dpi = 300)
