#-------------------
# Confluence Feeding
#-------------------

# confluence feeding overlaps with waters of cool SST, low Chl, 
# slow currents (relative to surrounding speeds)


rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(marmap)
}

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/KIPE/Ratmanoff, Kerguelen Islands chick-rearing tracks checked.rds")
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# isolate track
ind_02 <- trax %>% filter(individual_id == "02")

#-------------------------------------
# Plot over confluence zones
#-------------------------------------

# get bathymetry
e <- ext(trax)
bathy <- getNOAA.bathy(lon1 = e[1], lon2 = e[2], lat1 = e[3], lat2 = e[4], resolution = 1)
bathy <- as.xyz(bathy) %>%
  rename(lon = V1, lat = V2, depth = V3) %>%
  rast()

# get depth contours
contours <- as.contour(bathy, levels = seq(-4000, 0, 1000))
crs(contours) <- "epsg:4326"

# read in SST for 2020 tracks
sst <- rast("E:/Satellite_Data/daily/sst/sst_2020.nc")

# filter to best date
this_date <- as_date("2020-02-26")
sst <- sst[[time(sst) == this_date]]

# crop to track extent
sst <- crop(sst, ext(trax))

# create lines for this track
lines <- as.lines(ind_02)

# plot with each feeding event
p1 <- ggplot() + geom_spatraster(data = sst) +
  geom_spatvector(data = contours) +
  geom_spatvector(data = lines, lwd = 1.5, col = "grey90") + 
  geom_spatvector(data = ind_02 %>% filter(state == "ARS"),
                  col = "black", size = 2) +
  scale_fill_viridis_c(na.value = "white", limits = c(2.8, 7.2), name = "Sea Surface \nTemperature (°C)") +
  scale_color_manual(values = c("white", "grey80")) +
  scale_y_continuous(limits = c(-51.88, -49.04), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw()
p1

# read in SSH for 2020 tracks
ssh <- rast("E:/Satellite_Data/daily/ssh/ssh_2020.nc")

# filter to best date
this_date <- as_date("2020-02-26")
ssh <- ssh[[time(ssh) == this_date]]

# crop to track extent
ssh <- crop(ssh, ext(trax))

# plot with each feeding event
ggplot() + geom_spatraster(data = ssh) +
  geom_spatvector(data = contours) +
  geom_spatvector(data = lines, lwd = 1.5, col = "grey90") + 
  geom_spatvector(data = ind_02 %>% filter(state == "ARS"),
                  col = "black", size = 2) +
  scale_fill_viridis_c(na.value = "white") +
  scale_color_manual(values = c("white", "grey80")) +
  scale_y_continuous(limits = c(-51.88, -49.03), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw()

# read in uo and vo for 2020 tracks
uo <- rast("E:/Satellite_Data/daily/uo/uo_2020.nc")
vo <- rast("E:/Satellite_Data/daily/vo/vo_2020.nc")

# filter to best date
uo <- uo[[time(uo) == this_date]]
vo <- vo[[time(vo) == this_date]]

# crop to track extent
uo <- crop(uo, ext(trax))
vo <- crop(vo, ext(trax))

# calculate speed
speed <- sqrt(uo^2 + vo^2)

# plot with each feeding event
ggplot() + geom_spatraster(data = speed) +
  geom_spatvector(data = contours) +
  geom_spatvector(data = lines, lwd = 1.5, col = "grey90") + 
  geom_spatvector(data = ind_02 %>% filter(state == "ARS"),
                  col = "black", size = 2) +
  scale_fill_viridis_c(na.value = "white") +
  scale_color_manual(values = c("white", "grey80")) +
  scale_y_continuous(limits = c(-51.88, -49.03), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw()


# export temperature plot
ggsave("output/imagery/fig components/confluence KG.png", p1,
       width = 12, height = 6, dpi = 300)
