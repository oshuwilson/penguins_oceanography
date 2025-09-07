# KIPE Falklands

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)
library(terra)
library(tidyterra)
library(rnaturalearth)


# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/KIPE/Volunteer Beach, Falkland Islands late chick-rearing tracks checked.RDS")

# convert tracks to terra
tracks <- vect(tracks, geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")
# plot tracks
plot(tracks)

# get dates
unique(year(tracks$date))
unique(as_date(tracks$date))

# read in sst
sst <- rast("E:/Satellite_Data/daily/sst/sst_2008.nc")

# limit to early season
sst <- sst[[time(sst) < as_date("2008-02-01")]]

# crop to region of intensive activity
e <- ext(-62, -45, -60, -50)
sst <- crop(sst, e)

# # read in mld
# mld <- rast("E:/Satellite_Data/daily/mld/mld_2008.nc")
# mld <- mld[[time(mld) < as_date("2008-02-01")]]
# mld <- crop(mld, e)
# plot(mld[[9]])
# 
# # read in eddies
# eddies <- rast("E:/Satellite_Data/daily/eddies/eddies_2008.nc")
# eddies <- eddies[[time(eddies) < as_date("2008-02-01")]]
# eddies <- crop(eddies, e)
# plot(eddies[[15]])
# 
# # read in epipelagic nekton
# nekton <- rast("E:/Satellite_Data/daily/epipelagic_nekton/epipelagic_nekton_2008.nc")
# nekton <- nekton[[time(nekton) < as_date("2008-02-01")]]
# nekton <- crop(nekton, ext(-62, -45, -60, -53))
# plot(nekton[[10:15]])
# 
# 
# # read in chlorophyll
# chl <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2008_resampled.nc")
# chl <- chl[[time(chl) < as_date("2008-02-01")]]
# chl <- crop(chl, e)
# ggplot() +
#   geom_spatraster(data = chl[[1:10]]) +
#   scale_fill_viridis_c(trans = "log") +
#   facet_wrap(~lyr)

# plot a few SSTs
plot(sst[[11:15]])

# create lines from tracks
inds <- unique(tracks$individual_id)
trax <- tracks

for(ind in inds){
  
  ind_trax <- trax %>%
    filter(individual_id == ind)
  start_date <- min(ind_trax$date)
  ind_trax <- as.lines(ind_trax)
  ind_trax$date <- start_date
  ind_trax$individual_id <- ind
  if(ind == inds[1]){
    trax_lines <- ind_trax
  } else {
    trax_lines <- c(trax_lines, ind_trax)
  }
}
trax_lines <- vect(trax_lines)

# load in coastline of Falklands
fk <- ne_countries(scale = 10, returnclass = "sv") %>%
  crop(e)

# bounding box for eddy
bbox <- ext(-54, -50, -56.5, -54) %>%
  vect(crs = "epsg:4326")

# plot all tracks
all_tracks <- ggplot() + 
  geom_spatraster(data = sst[[10]]) + 
  geom_spatvector(data = fk) +
  geom_spatvector(data = trax_lines %>%
                    filter(!individual_id %in% c("Gus", "Leo", "Iona", "Caldera")), aes(col = individual_id), linewidth = 1) +
  geom_spatvector(data = bbox, fill = NA, col = "red4", linewidth = 2) +
  scale_x_continuous(limits = c(-61.96, -45.04), expand = c(0,0)) +
  scale_y_continuous(limits = c(-59.96, -50), expand = c(0,0)) +
  scale_fill_viridis_c(name = bquote(Sea ~ Surface ~ Temperature ~ (degree*C)),
                       na.value = "white") +
  scale_color_grey(name = "Individual ID", guide = "none") +
  theme_classic()
all_tracks

# plot Hansueli within bbox
hansueli <- ggplot() +
  geom_spatraster(data = sst[[10]]) + 
  geom_spatvector(data = fk) +
  geom_spatvector(data = trax_lines %>%
                    filter(individual_id == "Hansueli"), aes(col = individual_id), linewidth = 1) +
  scale_x_continuous(limits = c(-54.044, -50), expand = c(0,0)) +
  scale_y_continuous(limits = c(-56.543, -54), expand = c(0,0)) +
  scale_fill_viridis_c(name = bquote(Sea ~ Surface ~ Temperature ~ (degree*C)),
                       na.value = "white") +
  scale_color_grey(name = "Individual ID", guide = "none") +
  theme_classic()
hansueli

# save plots
ggsave("output/imagery/fig components/shoaling_all_tracks.png",
       all_tracks, width = 10, height = 8, dpi = 300, units = "in")
ggsave("output/imagery/fig components/shoaling_hansueli.png",
       hansueli, width = 10, height = 8, dpi = 300, units = "in")
time(sst[[10]])
