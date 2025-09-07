#----------------------------
# Create eddy trajectory plot
#----------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Predictor Data/eddy atlas")

{
  library(terra)
  library(tidyterra)
  library(sf)
  library(ncdf4)
  library(ggplot2)
  library(lubridate)
  library(marmap)
  library(rnaturalearth)
  library(rnaturalearthhires)
  library(rnaturalearthdata)
}

# read in tracks
tracks <- readRDS("~/OneDrive - University of Southampton/Documents/Chapter 02/output/hmm/hmm_tracks_by_colony/MAPE/Funk Beach, Marion Island incubation tracks checked.RDS")

# vectorise
trax <- vect(tracks, geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# get extent
e <- ext(25, 40, -53, -46)

# define area of interest
min_lat <- e[3] + 0.5
max_lat <- e[4] + 0.5
min_lon <- e[1] + 0.5
max_lon <- e[2] + 0.5

# correct longitude if below 0
if(min_lon < 0){
  min_lon <- 360 + min_lon
}
if(max_lon < 0){
  max_lon <- 360 + max_lon
}

# create minimum convex polygon
mcp <- convHull(trax)

# load in long cyclonic eddies
cyc_eddy <- nc_open("long_cyclones.nc")

# get center latitude, center longitude, radius, and time
lon <- as.data.frame(ncvar_get(cyc_eddy, "longitude"))
lat <- as.data.frame(ncvar_get(cyc_eddy, "latitude"))
rad <- as.data.frame(ncvar_get(cyc_eddy, "speed_radius"))
time <- as.data.frame(ncvar_get(cyc_eddy, "time"))
id <- as.data.frame(ncvar_get(cyc_eddy, "track"))

# join columns together
cyclones <- cbind(lon, lat, rad, time, id)
rm(list=setdiff(ls(), c("cyclones", "trax", "e")))
names(cyclones) <- c("lon", "lat", "radius", "time", "id")


# find cyclone IDs that cross area of interest
ids <- cyclones %>% filter(lat > min_lat & lat < max_lat & lon > min_lon  & lon < max_lon) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  select(id)

# only keep those cyclones
cyclones <- cyclones %>% filter(id %in% ids$id)
rm(ids)

# transform longitude
cyclones$lon <- ifelse(cyclones$lon > 180, cyclones$lon - 360, cyclones$lon)

# transform time
base_julian <- as_date("1950-01-01")
cyclones$time <- base_julian + days(cyclones$time)

# calculate relative time
cyclones <- cyclones %>%
  group_by(id) %>%
  mutate(rel_time = (time - min(time))/length(time))

# filter to those that took place in the date range of tracks
min_date <- as_date(min(trax$date))
max_date <- as_date(max(trax$date))
ids <- cyclones %>% filter(time >= min_date & time <= max_date) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  select(id)
cyclones <- cyclones %>% filter(id %in% ids$id)
rm(ids)

# plot
tracks <- trax %>% as.data.frame(geom = "XY")
cyclones %>% ggplot(aes(lon, lat)) +
  geom_point(data=tracks, aes(x, y), col="black") +
  geom_point(aes(col=as.factor(id))) +
  geom_path(aes(col=as.factor(id), group = id)) +
  coord_fixed() + 
  theme(legend.position = "none") 

# plot with relative time
cyclones %>% ggplot(aes(lon, lat)) +
  geom_point(aes(col=rel_time)) +
  geom_path(aes(col=rel_time, group = id)) +
  geom_point(data=tracks, aes(x, y), col="black") +
  coord_fixed() +
  scale_color_viridis_c()


#-------------------------
# Repeat for anticyclones
#-------------------------

# load in long anticyclonic eddies
anti_eddy <- nc_open("long_anticyclones.nc")

# get center latitude, center longitude, radius, and time
lon <- as.data.frame(ncvar_get(anti_eddy, "longitude"))
lat <- as.data.frame(ncvar_get(anti_eddy, "latitude"))
rad <- as.data.frame(ncvar_get(anti_eddy, "speed_radius"))
time <- as.data.frame(ncvar_get(anti_eddy, "time"))
id <- as.data.frame(ncvar_get(anti_eddy, "track"))

# join columns together
anticyclones <- cbind(lon, lat, rad, time, id)
rm(lon, lat, rad, time, id)
names(anticyclones) <- c("lon", "lat", "radius", "time", "id")

# find anticyclone IDs that cross area of interest
ids <- anticyclones %>% filter(lat > min_lat & lat < max_lat & lon > min_lon  & lon < max_lon) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  select(id)

# only keep those anticyclones
anticyclones <- anticyclones %>% filter(id %in% ids$id)
rm(ids)

# transform longitude
anticyclones$lon <- ifelse(anticyclones$lon > 180, anticyclones$lon - 360, anticyclones$lon)

# transform time
base_julian <- as_date("1950-01-01")
anticyclones$time <- base_julian + days(anticyclones$time)

# calculate relative time
anticyclones <- anticyclones %>%
  group_by(id) %>%
  mutate(rel_time = (time - min(time))/length(time))

# filter to those that took place in the date range of tracks
min_date <- as_date(min(trax$date))
max_date <- as_date(max(trax$date))
ids <- anticyclones %>% filter(time >= min_date & time <= max_date) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  select(id)
anticyclones <- anticyclones %>% filter(id %in% ids$id)
rm(ids)

# plot
tracks <- trax %>% as.data.frame(geom = "XY")
anticyclones %>% ggplot(aes(lon, lat)) +
  geom_point(data = tracks, aes(x, y), col = "black") +
  geom_point(aes(col=as.factor(id))) +
  geom_path(aes(col=as.factor(id), group = id)) +
  coord_fixed() +
  theme(legend.position = "none")

# plot with relative time
anticyclones %>% ggplot(aes(lon, lat)) +
  geom_point(data=tracks, aes(x, y), col="black") +
  geom_point(aes(col=rel_time)) +
  geom_path(aes(col=rel_time, group = id)) +
  coord_fixed() +
  scale_color_viridis_c()


#------------------
# Map Trajectories
#-----------------

# join cyclones and anticyclones
cyclones$type <- "Cyclones"
anticyclones$type <- "Anticyclones"
eddies <- rbind(anticyclones)

# cleanup
rm(list = setdiff(ls(), c("eddies", "trax", "e")))

# convert to terra
eddies <- vect(eddies, geom = c("lon", "lat"), crs = "epsg:4326")
e <- ext(trax) + ext(eddies)

# read in bathymetry
bathy <- getNOAA.bathy(lon1 = e[1], lon2 = e[2], lat1 = e[3], lat2 = e[4], resolution = 1)
bathy <- as.xyz(bathy) %>%
  rename(lon = V1, lat = V2, depth = V3) %>%
  rast()

# plot
ggplot(eddies) + 
  geom_spatraster(data = bathy) +
  geom_spatvector(data = trax) +
  geom_spatvector(aes(col = as.numeric(rel_time))) +
  facet_wrap(~type) +
  scale_color_viridis_c(end = 0.9) +
  scale_fill_gradient(low = "grey20", high = "white", limits = c(min(values(bathy)), 0), na.value = "white") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) 

# limit to only longer eddies
longeddies <- eddies %>%
  #filter(type == "Antiyclones") %>%
  group_by(id) %>%
  filter(dplyr::n() > 50) %>%
  ungroup()
nrow(longeddies)

mcp <- convHull(trax)
mcp <- buffer(mcp, 10000)
plot(mcp)

# create line for each eddy
ids <- unique(longeddies$id)
rm(eddy_lines)
for(this_id in ids){
  
  # filter to this eddy track
  this_eddy <- eddies %>% filter(id == this_id)
  
  # if no eddy points fall within mcp, skip
  within <- crop(this_eddy, mcp)
  if(nrow(within) == 0){
    next
  }

  
  # create line
  this_eddy_line <- as.lines(this_eddy)
  
  # assign ID
  this_eddy_line$id <- this_id
  
  # extract start and end point
  this_eddy_start <- this_eddy[1]
  this_eddy_end <- this_eddy[nrow(this_eddy)]
  
  # join all together
  if(!exists("eddy_lines")){
    eddy_lines <- this_eddy_line
    eddy_starts <- this_eddy_start
    eddy_ends <- this_eddy_end
  } else {
    eddy_lines <- c(eddy_lines, this_eddy_line)
    eddy_starts <- bind_spat_rows(eddy_starts, this_eddy_start)
    eddy_ends <- bind_spat_rows(eddy_ends, this_eddy_end)
  }
}

# make lines one spatvector
eddy_lines <- vect(eddy_lines)

# classify depth
mm1 <- c(-7000, -3500, 0,
         -3500, -2000, 1,
         -2000, -1000, 2,
         -1000, 0, 3)
bathy2 <- classify(bathy, mm1)
bathy2 <- as.factor(bathy2)
bathy2 <- droplevels(bathy2, level = c(8, 9))
plot(bathy2)

# add funk beach location
meta <- readRDS("~/OneDrive - University of Southampton/Documents/Chapter 02/data/metadata.RDS")
fb <- meta %>%
  filter(deployment_site == "Funk Beach, Marion Island") %>%
  slice(1) %>%
  select(deployment_decimal_longitude, deployment_decimal_latitude)
fb <- vect(fb, geom = c("deployment_decimal_longitude", "deployment_decimal_latitude"), crs = "epsg:4326")

# create bounding box
ext(trax)
bbox <- ext(37, 46, -52.5, -46)
bbox <- as.polygons(bbox)
crs(bbox) <- "epsg:4326"

# plot
p1 <- ggplot(eddy_lines) + 
  geom_spatraster(data = bathy2) +
  geom_spatvector(aes(col = as.factor(id)), lwd = 1.2) +
  geom_spatvector(data = eddy_starts, shape = 19, size = 2) +
  geom_spatvector(data = eddy_ends, shape = 15, size = 2) +
  geom_spatvector(data = bbox, lwd = 0.6, fill = NA, col = "black") +
  scale_color_viridis_d() +
  scale_fill_manual(values = c("grey90", "grey50", "grey30", "grey10"), na.value = "white",
                    name = "Depth (m)", 
                    labels = c("-7000 to -3500", "-3500 to -2000", "-2000 to -1000", "-1000 to 0", ""),
                    na.translate = F) +
  theme_classic() +
  scale_x_continuous(expand = c(0,0), limits = c(25, 50.1)) +
  scale_y_continuous(expand = c(0,0), limits = c(-52.9, -44.07)) +
  geom_spatvector(data = fb, size = 5, col = "black", fill = "red3", shape = 21) +
  guides(color = "none")
p1

# make lines for each track
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

# get outline of Marion
mar <- ne_countries(scale = 10, returnclass = "sv")
mar <- crop(mar, bbox)

# plot tracks and eddies in bounding box limits
p2 <- ggplot(eddy_lines) + 
  geom_spatraster(data = bathy2) +
  #geom_spatvector(aes(col = as.factor(id)), lwd = 1) +
  # geom_spatvector(data = eddy_starts, shape = 19, size = 2) +
  # geom_spatvector(data = eddy_ends, shape = 15, size = 2) +
  geom_spatvector(data = trax_lines, lwd = 1, col = "black") +
  geom_spatvector(data = mar, fill = "white") +
  scale_color_viridis_d() +
  scale_fill_manual(values = c("grey90", "grey50", "grey30", "grey10"), na.value = "grey10",
                    name = "Depth (m)", 
                    breaks = c(0, 2, 4, 6),
                    labels = c("-7000 to -3500", "-3500 to -2000", "-2000 to -1000", "-1000 to 0")) +
  theme_classic() +
  scale_x_continuous(expand = c(0,0), limits = c(37, 46)) +
  scale_y_continuous(expand = c(0,0), limits = c(-52.48, -46)) +
  geom_spatvector(data = fb, size = 4, col = "black", fill = "red3", shape = 21) +
  guides(color = "none")
p2

p1
p2

ggsave("~/OneDrive - University of Southampton/Documents/Chapter 02/output/imagery/fig components/funk_beach_eddies.png",
       p1, width = 12, height = 6.5, dpi = 300)
ggsave("~/OneDrive - University of Southampton/Documents/Chapter 02/output/imagery/fig components/funk_beach_tracks.png",
       p2, width = 8, height = 6.5, dpi = 300)

# length of eddy tracks
track_lengths <- perim(eddy_lines)
max(track_lengths)
min(track_lengths)
mean(track_lengths)
hist(track_lengths)
