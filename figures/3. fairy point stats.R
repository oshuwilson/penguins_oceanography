#-----------------------------------------------
# Quantify eddy impacts on nekton at Fairy Point
#-----------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

library(tidyverse)
library(terra)
library(tidyterra)
library(mgcv)
library(gratia)
library(cowplot)

# # read in background data for chick-rearing Fairy Point Macs
# bg <- readRDS("output/background/MAPE/Fairy Point, Bird Island early chick-rearing background.RDS")
# 
# # convert to terra
# bg <- vect(bg,
#            geom = c("x", "y"),
#            crs = "EPSG:4326")
# 
# # get years of background data
# years <- bg %>%
#   mutate(year = year(date)) %>%
#   pull(year) %>%
#   unique() %>%
#   sort()
# 
# # get dates
# dates <- bg %>%
#   mutate(date = as.Date(date)) %>%
#   pull(date) %>%
#   unique() %>%
#   sort()
# 
# # load extraction functions
# source("code/functions/extraction_functions.R")
# 
# # extract epipelagic nekton to background date
# bg_nek <- dynamic_extract("epipelagic_nekton", bg)
# 
# # reclass eddies between -1 and 1 to 0
# bg_nek <- bg_nek %>%
#   mutate(eddies = ifelse(eddies < 1 & eddies > -1, 0, eddies))
# 
# # plot eddy type vs nekton concentrations beyond the shelf
# ggplot(bg_nek %>% filter(depth > 0), aes(x = eddies, y = epipelagic_nekton)) +
#   geom_point() +
#   geom_smooth()
# 
# # extract chlorophyll
# bg_chl <- dynamic_chlorophyll("chlorophyll", bg)
# 
# # reclass eddies between -1 and 1 to 0
# bg_chl <- bg_chl %>%
#   mutate(eddies = ifelse(eddies < 1 & eddies > -1, 0, eddies))
# 
# # plot eddy type vs chlorophyll concentrations beyond the shelf
# ggplot(bg_chl, aes(x = eddies, y = chlorophyll)) +
#   geom_point() +
#   geom_smooth()
# 
# 
# # try with tracks
# # read in tracks
# tracks <- readRDS("output/hmm/hmm_tracks_by_colony/MAPE/Fairy Point, Bird Island early chick-rearing tracks checked.RDS")
# 
# # convert to terra
# tracks <- vect(tracks,
#                 geom = c("x", "y"),
#                 crs = "EPSG:6932") %>%
#   project("EPSG:4326")
# 
# # extract epipelagic nekton to tracks
# tracks_nek <- dynamic_extract("epipelagic_nekton", tracks)
# 
# # reclass eddies between -1 and 1 to 0
# tracks_nek <- tracks_nek %>%
#   mutate(eddies = ifelse(eddies < 1 & eddies > -1, 0, eddies))
# 
# # plot eddy type vs nekton concentrations beyond the shelf
# ggplot(tracks_nek %>% filter(depth > 0), aes(x = eddies, y = epipelagic_nekton)) +
#   geom_smooth()
# ggplot(tracks_nek, aes(x = state, y = epipelagic_nekton)) + geom_boxplot()
# 
# 
# # extract chlorophyll
# tracks_chl <- dynamic_chlorophyll("chlorophyll", tracks)
# 
# # reclass eddies between -1 and 1 to 0
# tracks_chl <- tracks_chl %>%
#   mutate(eddies = ifelse(eddies < 1 & eddies > -1, 0, eddies))
# 
# # plot eddy type vs chlorophyll concentrations beyond the shelf
# ggplot(tracks_chl %>% filter(state == "ARS"), aes(x = eddies, y = chlorophyll)) +
#   geom_point() +
#   geom_smooth()
# ggplot(tracks_chl, aes(x = state, y = chlorophyll)) + geom_boxplot()
# 
# # try with temperature
# ggplot(bg_nek %>% filter(year(date) == 2004), aes(x = eddies, y = epipelagic_nekton)) + 
#   geom_point() +
#   geom_smooth()
# ggplot(bg_chl, aes(x = eddies, y = curr)) + geom_smooth()
# 
# e <- ext(-42, -30, -55, -52)
# 
# # read in chlorophyll for 2004
# chl <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2012_resampled.nc")
# chl <- chl[[month(time(chl)) == 1]]
# chl <- crop(chl, e)
# plot(chl[[10:18]])
# plot(tracks_nek, add = T, pch = ".")
# 
# ggplot() +
#   geom_spatraster(data = chl) +
#   scale_fill_viridis_c(trans = "log") +
#   facet_wrap(~lyr)
# 
# # read in eddies for 2004
# eddies <- rast("E:/Satellite_Data/daily/eddies/eddies_2004.nc")
# eddies <- eddies[[time(eddies) > as_date("2004-01-01") & time(eddies) < as_date("2004-02-01")]]
# eddies <- crop(eddies, e)
# plot(eddies[[26:29]])
# 
# # read in sst for 2004
# sst <- rast("E:/Satellite_Data/daily/sst/sst_2004.nc")
# sst <- sst[[time(sst) > as_date("2004-01-01") & time(sst) < as_date("2004-02-01")]]
# sst <- crop(sst, e)
# plot(sst[[13:24]])


# # load bathymetry
# library(marmap)
# bathy <- getNOAA.bathy(-42, -30, -55, -52, resolution = 1, keep = TRUE)
# bathy <- as.xyz(bathy) %>%
#   rename(lon = V1, lat = V2, depth = V3) %>%
#   rast()
# plot(bathy)
# 
# # create bounding box around Bird Island shelf
# bbox <- ext(-40, -38, -54, -53) %>%
#   vect(crs = "epsg:4326") 
# plot(bbox, add = T)
# 
# # sample background points within bbox
# bg <- spatSample(bbox, 10000)
# plot(bg, add = T, pch = ".")
# 
# # extract bathymetry to background points
# bg$depth <- extract(bathy, bg, ID = F)
# ggplot(bg, aes(x=depth)) + geom_histogram()
# 
# # repeat plots with depth ranges
# plot(bathy)
# plot(bg %>% filter(depth < -800 & depth > -3400),
#      add = T, pch = ".")
# 
# # filter background points to this depth
# bg <- bg %>% filter(depth < -800 & depth > -3400)
# 
# # add days between 1 and 31 randomly
# bg$day <- sample(1:31, nrow(bg), replace = TRUE)
# 
# # make month 1
# bg$month <- 1
# 
# # sample year from years
# bg_99 <- bg %>% mutate(year = 1999)
# bg_03 <- bg %>% mutate(year = 2003)
# bg_04 <- bg %>% mutate(year = 2004)
# bg_05 <- bg %>% mutate(year = 2005)
# bg_12 <- bg %>% mutate(year = 2012)
# bg <- bind_spat_rows(bg_99, bg_03, bg_04, bg_05, bg_12)
# 
# # create date column
# bg$date <- as.Date(paste(bg$year, bg$month, bg$day, sep = "-"))
# 
# # extract eddies to background points
# bg <- dynamic_extract("eddies", bg)
# 
# # extract chlorophyll to bakcground points
# bg <- dynamic_chlorophyll("chl", bg)
# 
# # extract epipelagic nekton to background points
# bg <- dynamic_extract("epipelagic_nekton", bg)
# 
# # export background values
# saveRDS(bg, "exploration/stats/fairy_point_shelf.RDS")

# # create eddie background shading values
# ed_vals <- data.frame(type = as.factor(c("anticyclonic core", "anticyclonic periphery", "not eddy",
#                                          "cyclonic periphery", "cyclonic core")),
#                       x = c(-2.5, -1.5, 0, 1.5, 2.5),
#                       xmin = c(-3, -2, -1, 1, 2),
#                       xmax = c(-2, -1, 1, 2, 3), 
#                       ymin = 0, ymax = 1)
# 
# ggplot(bg, aes(x = e2, y = log(chl))) +
#   geom_tile(data = ed_vals, alpha = 0.4, height = Inf, width = c(1, 1, 2, 1, 1),
#             aes(x = x, y = 0.5, fill = type)) +
#   geom_smooth(col = "black", fill = "grey40") +
#   scale_fill_manual(values = c("red4", "red3",
#                                "steelblue4", "steelblue3",
#                                "grey80"),
#                     breaks = c("anticyclonic core", "anticyclonic periphery",
#                                "cyclonic core", "cyclonic periphery", "not eddy"),
#                     labels = c("Anticyclonic Core", "Anticyclonic Periphery",
#                                "Cyclonic Core", "Cyclonic Periphery", "Not Eddy"),
#                     name = "Eddy Zone") +
#   scale_x_continuous(expand = c(0,0), limits = c(-3, 3)) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Log Chlorophyll Concentration (mg/m³)") +
#   xlab("Relative Eddy Distance")


# read in background values
bg <- readRDS("exploration/stats/fairy_point_shelf.RDS")

# reclass eddies between -1 and 1 to 0
bg <- bg %>%
  mutate(ed2 = ifelse(eddies < 1 & eddies > -1, 0, eddies))

# plot eddy type vs chlorophyll concentrations beyond the shelf
ggplot(bg, aes(x = ed2, y = log(chl))) +
  geom_smooth()

# fit a GAM to the data
m1 <- gam(log(chl) ~ s(ed2), data = bg)

# get smooths
sm <- smooth_estimates(m1, n = 1000) %>%
  add_confint()

# average odds ratios over .5 intervals
sm2 <- sm %>%
  mutate(ed2 = cut(ed2, breaks = c(-3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.125, 
                                   0.125, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3))) %>%
  group_by(ed2) %>%
  summarise(.estimate = mean(.estimate),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci))

# remove two categories with no associated eddy region
sm2 <- sm2 %>% 
  filter(!ed2 %in% c("(-1,-0.125]", "(0.125,1]"))

# get upper and lower limits from cuts
sm2 <- sm2 %>%
  mutate(x_tmp = str_sub(ed2, 2, -2)) %>%
  separate(x_tmp, c("min", "max"), sep = ",") %>%
  mutate_at(c("min", "max"), as.numeric)

# rename odds ratio
sm2 <- sm2 %>%
  rename(OR = .estimate)

# if lower and upper ci envelope 1, set odds ratio to 1
sm2 <- sm2 %>%
  mutate(OR = ifelse(.lower_ci < 1 & .upper_ci > 1, 1, OR))

# read in blank odds ratios if some missing
blanksm <- readRDS("output/imagery/colony odds ratios/sm blank.RDS")
sm2 <- bind_rows(sm2, blanksm) %>%
  distinct(ed2, .keep_all = T)

# create cyclonic dataframe
cyclones <- sm2 %>% 
  filter(min >= 1) %>%
  arrange(desc(min))

# create factor for plotting
cyclones <- cyclones %>%
  mutate(band = factor(1:nrow(cyclones)))

# truncate odds ratios above 3
cyclones <- cyclones %>%
  mutate(OR = ifelse(OR > 3, 3, OR))

# create anticyclonic dataframe
anticyclones <- sm2 %>% 
  filter(max <= -1) %>%
  arrange(min)

# create factor for plotting
anticyclones <- anticyclones %>%
  mutate(band = factor(1:nrow(anticyclones)))

# truncate odds ratios above 3
anticyclones <- anticyclones %>%
  mutate(OR = ifelse(OR > 3, 3, OR))

# get estimate for no eddy zone
no_eddy <- sm2 %>%
  filter(ed2 == "(-0.125,0.125]") %>%
  select(OR) %>%
  pull()

# get min and max values
min_val <- min(sm2$OR, na.rm = TRUE)
max_val <- max(sm2$OR, na.rm = TRUE)

# take biggest of both absolute values
abs_val <- max(abs(min_val), abs(max_val))

# change min and max to absolute
min_val <- -abs_val
max_val <- abs_val

# delta for plotting
d <- "Delta"

# plot cyclones
p2 <- ggplot(cyclones, aes(x = band, fill = as.numeric(OR))) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 8.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = no_eddy,
                       limits = c(min_val, max_val),
                       na.value = "grey") +
  ggtitle("Cyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5)) +
  labs(fill = bquote(.(sym(d)) ~ Log ~ Chlorophyll ~ Concentration ~ (mg ~ m^-3))) +
  theme(legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top")) 
p2

# plot anticyclones
p3 <- ggplot(anticyclones, aes(x = band, fill = as.numeric(OR))) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 8.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = no_eddy,
                       limits = c(min_val, max_val),
                       na.value = "grey") +
  ggtitle("Anticyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5)) +
  labs(fill = bquote(.(sym(d)) ~ Log ~ Chlorophyll ~ Concentration ~ (mg ~ m^-3))) +
  theme(legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top"))
p3

# export p2 and p3
ggsave("output/imagery/fig components/fairy_point_cyclones.png", p2, 
       width = 6, height = 6, dpi = 300)
ggsave("output/imagery/fig components/fairy_point_anticyclones.png", p3,
       width = 6, height = 6, dpi = 300)
