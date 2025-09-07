#---------------------------------------------------
# Explore Eddies around Pointe Geologie
#---------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)
library(mgcv)
library(gratia)

# read in tracks
tracks <- readRDS("output/auger_extractions/ADPE/Pointe Geologie_incubation_extracted.RDS")

# # convert to terra
# tracks <- vect(tracks,
#                geom = c("x", "y"),
#                crs = "EPSG:4326")

# limit to 2017
tracks <- tracks %>% 
  filter(year(date) == 2017)

# # load extraction functions
# source("code/functions/extraction_functions.R")
# 
# # extract chlorophyll
# tracks <- dynamic_chlorophyll("chlorophyll", tracks)
# 
# # plot 2017 chlorophyll
# chl <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2017_resampled.nc")
# chl <- chl[[time(chl) == as_date("2017-12-15")]]
# chl <- crop(chl, ext(tracks))
# 
# ggplot(chl) + 
#   geom_spatraster(data = chl) +
#   scale_fill_viridis_c(trans = "log")
# 
# # extract epipelagic_nekton
# tracks <- dynamic_extract("epipelagic_nekton", tracks)
# 
# # plot 2017 epipelagic nekton
# nekton <- rast("E:/Satellite_Data/daily/epipelagic_nekton/epipelagic_nekton_2017.nc")
# nekton <- nekton[[time(nekton) == as_date("2017-12-13")]]
# nekton <- crop(nekton, ext(tracks))
# 
# # plot 2017 upper_mesopelagic_nekton
# nekton <- rast("E:/Satellite_Data/daily/upper_meso_nekton/upper_meso_nekton_2017.nc")
# nekton <- nekton[[time(nekton) == as_date("2017-12-13")]]
# nekton <- crop(nekton, ext(tracks))
# 
# # plot 2017 upper_meso_mig_nekton
# nekton <- rast("E:/Satellite_Data/daily/upper_mig_meso_nekton/upper_mig_meso_nekton_2017.nc")
# nekton <- nekton[[time(nekton) == as_date("2017-12-13")]]
# nekton <- crop(nekton, ext(tracks))
# 
# ggplot() +
#   geom_spatraster(data = nekton) +
#   scale_fill_viridis_c(trans = "sqrt")
# 
# 
# # plot 2017 salinity
# sal <- rast("E:/Satellite_Data/daily/sal/sal_2017.nc")
# sal <- sal[[time(sal) == as_date("2017-12-30")]]
# sal <- crop(sal, ext(137, 148, -67, -62))
# ggplot() +
#   geom_spatraster(data = sal) +
#   scale_fill_viridis_c()
# 
# # plot 2017 ssh
# ssh <- rast("E:/Satellite_Data/daily/ssh/ssh_2017.nc")
# ssh <- ssh[[time(ssh) == as_date("2017-12-31")]]
# ssh <- crop(ssh, ext(137, 148, -67, -62))
# ggplot() +
#   geom_spatraster(data = ssh) +
#   scale_fill_viridis_c()
# 
# # plot 2017 eddies
# eddies <- rast("E:/Satellite_Data/daily/eddies_auger/eddies_auger_2017.nc")
# eddies <- eddies[[time(eddies) >= as_date("2017-12-13")]]
# eddies <- crop(eddies, ext(137, 148, -67, -62))
# plot(eddies)
# ggplot() +
#   geom_spatraster(data = eddies) +
#   scale_fill_viridis_c() +
#   facet_wrap(~lyr)

# reclass all eddies between -1 and 1 to 0
tracks <- tracks %>% 
  mutate(eddies_auger = ifelse(eddies_auger <= -1 | eddies_auger >= 1, eddies_auger, 0))

# create eddie background shading values
ed_vals <- data.frame(type = as.factor(c("anticyclonic core", "anticyclonic periphery", "not eddy",
                                         "cyclonic periphery", "cyclonic core")),
                      x = c(-2.5, -1.5, 0, 1.5, 2.5),
                      xmin = c(-3, -2, -1, 1, 2),
                      xmax = c(-2, -1, 1, 2, 3), 
                      ymin = 0, ymax = 1)

# plot salinity vs eddies
p1 <- ggplot(tracks, aes(x = eddies_auger, y = sal)) + 
  geom_tile(data = ed_vals, alpha = 0.4, height = Inf, width = c(1, 1, 2, 1, 1),
            aes(x = x, y = 0.5, fill = type)) +
  geom_smooth(col = "black", fill = "grey40") +
  scale_fill_manual(values = c("red4", "red3",
                               "steelblue4", "steelblue3",
                               "grey80"),
                    breaks = c("anticyclonic core", "anticyclonic periphery",
                               "cyclonic core", "cyclonic periphery", "not eddy"),
                    labels = c("Anticyclonic Core", "Anticyclonic Periphery",
                               "Cyclonic Core", "Cyclonic Periphery", "Not Eddy"),
                    name = "Eddy Zone") +
  coord_cartesian(ylim = c(33.65, 33.95)) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Salinity (PSU)") +
  xlab("Relative Eddy Distance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p1

# rename eddies_auger to ed2
tracks <- tracks %>%
  rename(ed2 = eddies_auger)

# plot eddy type vs chlorophyll concentrations beyond the shelf
ggplot(tracks, aes(x = ed2, y = sal)) +
  geom_smooth()

# fit a GAM to the data
m1 <- gam(sal ~ s(ed2), data = tracks)

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

# subtract no_eddy stat from all odds ratios
sm2 <- sm2 %>%
  mutate(OR = OR - no_eddy)

# get min and max values
min_val <- min(sm2$OR, na.rm = TRUE)
max_val <- max(sm2$OR, na.rm = TRUE)

# take biggest of both absolute values
abs_val <- max(abs(min_val), abs(max_val))

# change min and max to absolute
min_val <- -abs_val
max_val <- abs_val

# plot cyclones
p2 <- ggplot(cyclones, aes(x = band, fill = as.numeric(OR))) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 8.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 0,
                       limits = c(min_val, max_val),
                       na.value = "grey") +
  ggtitle("Cyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5)) +
  labs(fill = expression(Delta~"Salinity (PSU)")) +
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
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 0,
                       limits = c(min_val, max_val),
                       na.value = "grey") +
  ggtitle("Anticyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5)) +
  labs(fill = expression(Delta~"Salinity (PSU)")) +
  theme(legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top"))
p3 <- p3 + theme(legend.key.size = unit(1, "cm"))

# export plots
#ggsave("output/imagery/fig components/pointe_geologie_salinity.png", p1, width = 10, height = 7, dpi = 300)
ggsave("output/imagery/fig components/pointe_geologie_cyclones.png", p2, width = 6, height = 6, dpi = 300)
ggsave("output/imagery/fig components/pointe_geologie_anticyclones.png", p3, width = 6, height = 6, dpi = 300)
