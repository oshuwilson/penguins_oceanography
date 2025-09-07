#---------------------------------------------------
# Explore Eddies around Funk Beach
#---------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)
library(terra)
library(tidyterra)

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/MAPE/Funk Beach, Marion Island incubation tracks checked.RDS")

# # convert to terra
tracks <- vect(tracks,
               geom = c("x", "y"),
               crs = "EPSG:6932") %>%
  project("epsg:4326")

# find incubation years
years <- unique(year(tracks$date))

# # read in nekton for 2018
# nek <- rast("E:/Satellite_Data/daily/upper_mig_meso_nekton/upper_mig_meso_nekton_2018.nc")
# nek <- nek[[time(nek) >= as_date("2018-12-13")]]
# nek <- crop(nek, ext(tracks))
# plot(nek[[10]])
# 
# # read in epipelagic_nekton for 2018
# epinekt <- rast("E:/Satellite_Data/daily/epipelagic_nekton/epipelagic_nekton_2018.nc")
# epinekt <- epinekt[[time(epinekt) >= as_date("2018-12-13")]]
# epinekt <- crop(epinekt, ext(tracks))
# plot(epinekt[[10]])
# 
# # read in chlorophyll for 2018
# chl <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2018_resampled.nc")
# chl <- chl[[time(chl) >= as_date("2018-12-13")]]
# chl <- crop(chl, ext(bbox))
# plot(chl[[1:9]])
# 
# # read in eddies for 2018
# eddies <- rast("E:/Satellite_Data/daily/eddies/eddies_2018.nc")
# eddies <- eddies[[time(eddies) >= as_date("2018-12-13")]]
# eddies <- crop(eddies, ext(tracks))
# plot(eddies)

# # create channel feature
# chan <- points %>%
#   as.data.frame(geom = "XY") %>%
#   mutate(channel = ifelse(depth < -3500 & y < -46, "channel", "outside"))
# ggplot(chan, aes(x=channel, y=log(chlorophyll))) + geom_boxplot(outliers = F)
# 
# # read in eddies for each year
# eddies18 <- rast("E:/Satellite_Data/daily/eddies/eddies_2018.nc") 
# eddies18 <- eddies18[[time(eddies18) %in% unique(as_date(tracks$date))]]
# eddies18 <- crop(eddies18, ext(bbox))
# plot(eddies18)
# 
# eddies19 <- rast("E:/Satellite_Data/daily/eddies/eddies_2019.nc")
# eddies19 <- eddies19[[time(eddies19) %in% unique(as_date(tracks$date))]]
# eddies19 <- crop(eddies19, ext(bbox))
# plot(eddies19)
# 
# eddies11 <- rast("E:/Satellite_Data/daily/eddies/eddies_2011.nc")
# eddies11 <- eddies11[[time(eddies11) %in% unique(as_date(tracks$date))]]
# eddies11 <- crop(eddies11, ext(bbox))
# plot(eddies11)
# 
# eddies12 <- rast("E:/Satellite_Data/daily/eddies/eddies_2012.nc")
# eddies12 <- eddies12[[time(eddies12) %in% unique(as_date(tracks$date))]]
# eddies12 <- crop(eddies12, ext(bbox))
# plot(eddies12)
# 
# 
# # combine eddies11 and eddies12
# eddies11_12 <- c(eddies11, eddies12)
# 
# # split eddies19 by month and combine january dates with eddies18
# eddies19_jan <- eddies19[[month(time(eddies19)) == 1]]
# eddies19_20 <- eddies19[[month(time(eddies19)) != 1]]
# eddies18_19 <- c(eddies18, eddies19_jan)
# 
# # revalue eddies to 1 or 0
# m1 <- matrix(c(-3, -0.01, 1,
#                0, 0, 0,
#                0.01, 3, 1),
#              nrow = 3,
#              byrow = T)
# 
# eddies11_12 <- classify(eddies11_12, m1, include.lowest = TRUE, right = TRUE)
# eddies18_19 <- classify(eddies18_19, m1, include.lowest = TRUE, right = TRUE)
# eddies19_20 <- classify(eddies19_20, m1, include.lowest = TRUE, right = TRUE)
# 
# # calculate eddy prevalence
# prev11_12 <- sum(eddies11_12)/nlyr(eddies11_12)
# prev18_19 <- sum(eddies18_19)/nlyr(eddies18_19)
# prev19_20 <- sum(eddies19_20)/nlyr(eddies19_20)
# plot(prev11_12)
# 
# # read in chlorophyll concentrations over the same timeframe
# chl11 <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2011_resampled.nc")
# chl12 <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2012_resampled.nc")
# chl18 <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2018_resampled.nc")
# chl19 <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2019_resampled.nc")
# 
# # subset chlorophyll to the same dates as eddies
# chl11 <- chl11[[time(chl11) %in% unique(as_date(tracks$date))]]
# chl12 <- chl12[[time(chl12) %in% unique(as_date(tracks$date))]]
# chl18 <- chl18[[time(chl18) %in% unique(as_date(tracks$date))]]
# chl19 <- chl19[[time(chl19) %in% unique(as_date(tracks$date))]]
# 
# # crop chlorophyll to the bounding box
# chl11 <- crop(chl11, ext(bbox))
# chl12 <- crop(chl12, ext(bbox))
# chl18 <- crop(chl18, ext(bbox))
# chl19 <- crop(chl19, ext(bbox))
# 
# 
# # combine chlorophyll layers
# chl11_12 <- c(chl11, chl12)
# 
# # split chl19 by month and combine january dates with chl18
# chl19_jan <- chl19[[month(time(chl19)) == 1]]
# chl19_20 <- chl19[[month(time(chl19)) != 1]]
# chl18_19 <- c(chl18, chl19_jan)
# 
# # calculate average chlorophyll
# chl11_12_avg <- sum(chl11_12) / nlyr(chl11_12)
# plot(chl11_12_avg)
# plot(prev11_12)
# 
# # calculate average chlorophyll for 2018 and 2019
# chl18_19_avg <- sum(chl18_19) / nlyr(chl18_19)
# ggplot() +
#   geom_spatraster(data = chl18_19_avg) +
#   scale_fill_viridis_c(trans = "log")
# plot(prev18_19)
# 
# # calculate average chlorophyll for 2019 and 2020
# chl19_20_avg <- sum(chl19_20) / nlyr(chl19_20)
# ggplot() +
#   geom_spatraster(data = chl19_20_avg) +
#   scale_fill_viridis_c(trans = "log")
# plot(prev19_20)
# 
# # sample prev18_19 and chl18_19_avg
# points$prev18_19 <- extract(prev18_19, points, ID = F)
# points$chl18_19_avg <- extract(chl18_19_avg, points, ID = F)
# 
# ggplot(points %>% filter(depth < -3500), aes(x = prev18_19, y = log(chl18_19_avg))) + 
#   geom_point() +
#   geom_smooth()
# 
# 
# # sample prev19_20 and chl19_20_avg
# points$prev19_20 <- extract(prev19_20, points, ID = F)
# points$chl19_20_avg <- extract(chl19_20_avg, points, ID = F)
# 
# ggplot(points %>% filter(depth < -3500), aes(x = prev19_20, y = log(chl19_20_avg))) + 
#   geom_point() +
#   geom_smooth()
# 
# 
# # sample prev11_12 and chl11_12_avg
# points$prev11_12 <- extract(prev11_12, points, ID = F)
# points$chl11_12_avg <- extract(chl11_12_avg, points, ID = F)
# 
# ggplot(points %>% filter(depth < 0), aes(x = prev11_12, y = log(chl11_12_avg))) + 
#   geom_point() +
#   geom_smooth()
# 
# # combine all prevalences
# prevs <- points %>%
#   as.data.frame() %>%
#   pivot_longer(cols = starts_with("prev"), 
#                names_to = "year", 
#                values_to = "eddy_prevalence") %>%
#   pivot_longer(cols = starts_with("chl1"), 
#                names_to = "year_chl", 
#                values_to = "chl_avg")
# 
# ggplot(prevs %>% filter(depth < -3500), aes(x = eddy_prevalence, y = log(chl_avg))) + 
#   geom_smooth()

# create bounding box to include tracks + ocean rises
bbox <- ext(37, 46, -53, -44) %>%
   vect(crs = "EPSG:4326")

# # load bathymetry
library(marmap)
bathy <- getNOAA.bathy(ext(bbox)[1], ext(bbox)[2], 
                       ext(bbox)[3], ext(bbox)[4], 
                       resolution = 1)
bathy <- as.xyz(bathy) %>%
  rename(lon = V1, lat = V2, depth = V3) %>%
  rast()
plot(bathy)

# # sample 10000 points from the bounding box
# points <- spatSample(bbox, 50000)
# 
# # extract depth to points
# points$depth <- extract(bathy, points, ID = F)
# ggplot(points, aes(x=depth)) + geom_histogram()
# 
# # load extraction functions
# source("code/functions/extraction_functions.R")
# 
# # append time to points from tracks
# points$date <- sample(tracks$date, 
#                       nrow(points), 
#                       replace = T)
# 
# # extract chlorophyll
# points <- dynamic_chlorophyll("chlorophyll", points)
# 
# # extract eddies
# points <- dynamic_extract("eddies", points)
# 
# # revalue eddies 
# points <- points %>%
#   mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))
# 
# # export points
#saveRDS(points, "exploration/stats/funk_beach.rds")
# 
# # extract nekton
# points <- dynamic_extract("epipelagic_nekton", points)
# points <- dynamic_extract("upper_mig_meso_nekton", points)
# points <- dynamic_extract("lower_hmig_meso_nekton", points)

# read in points
points <- readRDS("exploration/stats/funk_beach.rds")

plot(bathy)
plot(points %>% filter(depth < -3800), pch = ".", add = T, col = "white")  

ggplot(points %>% as.data.frame(geom= "XY") %>% filter(depth < 0 & y < -46 & y > -51.5), 
       aes(x=ed2, y=epipelagic_nekton)) +
  geom_smooth()


# plot channel vs productivity
ggplot(points %>% filter(depth < 0), aes(x=depth, y = log(chlorophyll))) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5))

# create channel background shading values
chan_vals <- data.frame(type = as.factor(c("deep channel", "ocean rise")),
                      x = c(-4500, - 1750),
                      xmin = c(-5500, -3500),
                      xmax = c(-3500, 0), 
                      ymin = -1.7, ymax = -1)

# plot channel vs productivity
p1 <- ggplot(points %>% filter(depth < 0), aes(x=depth, y = log(chlorophyll))) + 
  geom_tile(data = chan_vals, alpha = 0.4, height = Inf, width = c(2000, 3500),
            aes(x = x, y = -1.5, fill = type)) +
  geom_smooth(col = "black", fill = "grey40",
              method = "gam", formula = y ~ s(x, k = 5)) +
  scale_fill_manual(values = c("steelblue4",
                               "steelblue1"),
                    labels = c("Deep Channel",
                               "Ocean Rise"),
                    name = "Depth Category") +
  coord_cartesian(xlim = c(min(points$depth), 0)) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  ylab(bquote(Log ~ Chlorophyll ~ Concentration ~ (mg ~ m^-3))) +
  xlab("Seafloor Depth (m)")
p1

# export plot
ggsave("output/imagery/fig components/funk_beach_chlorophyll_depth.png",
       p1, width = 8, height = 6, dpi = 300)
