#-------------------------
# Interannual Variability
#-------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(marmap)
  library(rnaturalearth)
}

# 1. Incubating Macs from Fairy Point

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/MAPE/Fairy Point, Bird Island incubation tracks checked.rds")

# calculate max distance to colony for each individual
dists <- tracks %>%
  group_by(individual_id) %>%
  summarise(max_dist = max(dist_to_colony))

# calculate median distance of eddy foraging events
eddy_dists <- tracks %>%
  filter(eddies <= -1 | eddies >= 1) %>%
  filter(state == "ARS") %>%
  group_by(individual_id) %>%
  summarise(eddy_dist = median(dist_to_colony))

# calculate earliest eddy ARS encounter time
tracks <- tracks %>%
  group_by(individual_id) %>%
  mutate(init = if_else(state == "ARS" & eddies <= -1 | state == "ARS" & eddies >= 1, 1, 0)) %>%
  mutate(portion = cumsum(init))

# when portion reaches 10 hours of eddy-centric foraging, get time since start
eddy_time <- tracks %>%
  group_by(individual_id) %>%
  filter(portion == 20) %>%
  slice(1) %>%
  select(individual_id, time_since_start)

# join together
all <- dists %>%
  left_join(eddy_dists) %>%
  left_join(eddy_time)

# linear model
m1 <- lm(max_dist ~ time_since_start, data = all)
summary(m1)

# plot
ggplot(all, aes(x = time_since_start, y = max_dist)) +
  geom_point(col = "darkred") +
  geom_smooth(method = "lm", col = "black") +
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) 

# convert to terra
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# eddy foraging events
trax$init <- as.factor(trax$init)
eddy_events <- trax %>% filter(init == 1)

# read in fronts
fronts <- sf::read_sf("exploration/colony summaries/south georgia/sg_polar_fronts.shp")
fronts <- vect(fronts)
fronts <- project(fronts, "epsg:4326")
e <- ext(-46, -35, -55.5, -48)
fronts <- crop(fronts, e)

# load bathymetry
bathy <- getNOAA.bathy(lon1 = e[1], lon2 = e[2], lat1 = e[3], lat2 = e[4], resolution = 1)
bathy <- as.xyz(bathy) %>%
  rename(lon = V1, lat = V2, depth = V3) %>%
  rast()
plot(bathy)

# load coastline
coast <- ne_countries(scale = 10, returnclass = "sv")
coast <- crop(coast, e)

# create lines for each individual
inds <- unique(trax$individual_id)
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

# change front names
fronts <- fronts %>%
  filter(name != "Southern Antarctic Circumpolar Current Front (sACCf)") %>%
  mutate(name = ifelse(name == "Polar Front (PF)", "Polar Front", "Subantarctic Front"))

# Fairy Point colony location
meta <- readRDS("data/metadata.RDS")
fp <- meta %>% 
  filter(deployment_site == "Fairy Point, Bird Island") %>%
  slice(1) %>%
  select(deployment_decimal_longitude, deployment_decimal_latitude)
fp <- vect(fp, geom = c("deployment_decimal_longitude", "deployment_decimal_latitude"), crs = "epsg:4326")

# plot
p1 <- ggplot() +
  geom_spatraster(data = bathy) +
  geom_spatvector(data = fronts, lwd = 1.5, aes(linetype = name)) +
  geom_spatvector(data = coast, fill = "white") +
  geom_spatvector(data = trax_lines, col = "grey40", lwd = 1) +
  geom_spatvector(data = eddy_events, col = "goldenrod") +
  geom_spatvector(data = fp, size = 3, col = "black", fill = "red3", shape = 21) +
  scale_fill_gradient(limits = c(-6600, 0), name = "Depth (m)", 
                      high = "#eaeff3", low = "steelblue4", na.value = "#eaeff3") +
  facet_wrap(~year(date)) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_linetype(name = "Fronts")
p1

ggsave("output/imagery/fig components/interannual_fairy_point.png", p1, width = 10, height = 8, units = "in")

# calculate trip length for each individuals
lengths <- data.frame(track_len = perim(trax_lines), 
                      individual_id = trax_lines$individual_id)

# get year for each individual
ind_years <- tracks %>%
  group_by(individual_id) %>%
  summarise(year = min(year(date)))

# merge
lengths <- lengths %>%
  left_join(ind_years) %>%
  mutate(year = as.factor(year))

# plot
lengths %>%
  ggplot(aes(x = year, y = track_len)) +
  geom_boxplot()

# anova
m1 <- aov(track_len ~ year, data = lengths)
summary(m1)
plot(m1)

# group all years except 2000
lengths <- lengths %>%
  mutate(yeardiff = ifelse(year == 2000, "2000", "Other"))

# anova 2
m2 <- aov(track_len ~ yeardiff, data = lengths)
summary(m2)

# prep data for wilcox and t test
data2000 <- lengths %>%
  filter(year == 2000) %>%
  pull(track_len)
data_other <- lengths %>%
  filter(year != 2000) %>%
  pull(track_len)

# run t-test and wilcox test
t.test(data2000, data_other)
wilcox.test(data2000, data_other)


# mean length for each year
lengths %>%
  group_by(year) %>%
  summarise(mean_len = mean(track_len)/1000,
            min_len = min(track_len)/1000,
            max_len = max(track_len/1000))

# for 2000 vs other years
lengths %>%
  group_by(yeardiff) %>%
  summarise(mean_len = mean(track_len)/1000,
            min_len = min(track_len)/1000,
            max_len = max(track_len/1000),
            lq = quantile(track_len, 0.25)/1000,
            uq = quantile(track_len, 0.75)/1000)


# calculate track duration for each individual
durations <- tracks %>%
  group_by(individual_id) %>%
  filter(!is.na(date)) %>%
  summarise(max_date = max(date),
            min_date = min(date)) %>%
  mutate(duration = difftime(max_date, min_date, units = "days"),
         year = year(min_date))


# plot track duration
durations %>%
  ggplot(aes(x = as.factor(year), y = duration)) +
  geom_boxplot() +
  labs(x = "Year", y = "Duration (days)") +
  theme_bw()

# calculate mean and IQR for 2000 vs other years
durations <- durations %>%
  mutate(yeardiff = ifelse(year == 2000, "2000", "Other"))

durations %>%
  group_by(yeardiff) %>%
  summarise(mean_duration = mean(duration),
            min_duration = min(duration),
            max_duration = max(duration),
            lq = quantile(duration, 0.25),
            uq = quantile(duration, 0.75))

# prep data for t-test and wilcox test
durations <- durations %>%
  mutate(duration = as.numeric(duration))
data2000 <- durations %>%
  filter(year == 2000) %>%
  pull(duration)
data_other <- durations %>%
  filter(year != 2000) %>%
  pull(duration)

t.test(data2000, data_other)
wilcox.test(data2000, data_other)


# plot over eddy activity for each year

# list dates
dates <- unique(as_date(tracks$date))

# list years
years <- c(1999, 2000, 2001, 2011)

# for each year
year <- 2011

# read in eddies for this year
eddies <- rast(paste0("E:/Satellite_Data/daily/eddies/eddies_", year, ".nc"))

# limit to dates of tracks
eddies <- eddies[[time(eddies) %in% dates]]

# crop to extent
eddies <- crop(eddies, e)

# classify eddies as eddy (below -1/above 1) or non-eddy
m1 <- matrix(c(-3, -1, 1,
             -1, 1, 0,
             1, 3, 1),
             nrow = 3,
             byrow = T)

eds <- classify(eddies, m1, include.lowest = TRUE, right = TRUE)

# calculate sum of eds
eds_sum <- sum(eds)/nlyr(eds)
plot(eds_sum)

# plot the tracks for this year
tracks_year <- trax_lines %>%
  filter(year(date) == year)
plot(tracks_year, add = T)

# bank this year
#eds_1999 <- eds_sum
#eds_2000 <- eds_sum
#eds_2001 <- eds_sum
#eds_2011 <- eds_sum

# plot each eddy year one by one
p1999 <- ggplot() +
  geom_spatraster(data = eds_1999) +
  geom_spatvector(data = coast, fill = "white") +
  geom_spatvector(data = trax_lines %>% filter(year(date) == 1999), col = "grey75", lwd = 1) +
  geom_spatvector(data = eddy_events %>% filter(year(date) == 1999), col = "goldenrod") +
  geom_spatvector(data = fp, size = 3, col = "black", fill = "red3", shape = 21) +
  scale_fill_gradient(name = "Proportion of days\nwith eddy presence") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
p1999

p2000 <- ggplot() +
  geom_spatraster(data = eds_2000) +
  geom_spatvector(data = coast, fill = "white") +
  geom_spatvector(data = trax_lines %>% filter(year(date) == 2000), col = "grey75", lwd = 1) +
  geom_spatvector(data = eddy_events %>% filter(year(date) == 2000), col = "goldenrod") +
  geom_spatvector(data = fp, size = 3, col = "black", fill = "red3", shape = 21) +
  scale_fill_gradient(name = "Proportion of days\nwith eddy presence") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
p2000

p2001 <- ggplot() +
  geom_spatraster(data = eds_2001) +
  geom_spatvector(data = coast, fill = "white") +
  geom_spatvector(data = trax_lines %>% filter(year(date) == 2001), col = "grey75", lwd = 1) +
  geom_spatvector(data = eddy_events %>% filter(year(date) == 2001), col = "goldenrod") +
  geom_spatvector(data = fp, size = 3, col = "black", fill = "red3", shape = 21) +
  scale_fill_gradient(name = "Proportion of days\nwith eddy presence") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
p2001

p2011 <- ggplot() +
  geom_spatraster(data = eds_2011) +
  geom_spatvector(data = coast, fill = "white") +
  geom_spatvector(data = trax_lines %>% filter(year(date) == 2011), col = "grey75", lwd = 1) +
  geom_spatvector(data = eddy_events %>% filter(year(date) == 2011), col = "goldenrod") +
  geom_spatvector(data = fp, size = 3, col = "black", fill = "red3", shape = 21) +
  scale_fill_gradient(name = "Proportion of days\nwith eddy presence") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
p2011

# combine plots
library(cowplot)

# get legend
legend <- get_legend(p1999)

# remove legends from plots
p1999 <- p1999 + theme(legend.position = "none")
p2000 <- p2000 + theme(legend.position = "none")
p2001 <- p2001 + theme(legend.position = "none")
p2011 <- p2011 + theme(legend.position = "none")

# combine
plots <- plot_grid(p1999, p2000, p2001, p2011, ncol = 2, nrow = 2, align = "hv", 
                   labels = c("1999", "2000", "2001", "2011"),
                   label_x = 0.1, label_y = 0.97, label_colour = "white")
plots

# add legend
final_plot <- plot_grid(plots, legend, ncol = 2, rel_widths = c(1, 0.16))
final_plot

ggsave("output/imagery/fig components/8. interannual_fairy_point_eddies.png", final_plot, width = 9.5, height = 8, dpi = 300)




