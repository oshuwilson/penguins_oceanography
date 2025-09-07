#---------------------
# Thermocline Shoaling
#---------------------

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
tracks <- tracks %>%
  filter(eddies <= -1 | eddies >= 1)
ars <- tracks %>%
  filter(state == "ARS")

# plot mld against eddy zone
ggplot(tracks, aes(x = eddies)) + geom_histogram()
ggplot(trax, aes(x = year(date), fill = individual_id)) + geom_histogram()

# read in nekton
nek <- rast("E:/Satellite_Data/daily/lower_hmig_meso_nekton/lower_hmig_meso_nekton_2013.nc")

# limit to dates of trax
nek <- nek[[time(nek) %in% unique(as_date(trax$date))]]

# crop to track extent
nek <- crop(nek, ext(trax))
nek_avg <- mean(nek, na.rm = TRUE)

# plot
ggplot() + geom_spatraster(data = nek_avg) + 
  geom_spatvector(data = trax, size = 0.5)

# plot for each day
dates <- unique(as_date(trax$date))
dates <- dates[year(dates) == 2013]

ext(trax)
ggplot(trax) + geom_spatvector(aes(col = individual_id))
sub_trax <- trax %>%
  filter(year(date) == 2013)
ggplot(sub_trax, aes(x=date)) + geom_histogram(aes(fill = individual_id))
plot(sub_trax)

# read in temperature at depths
temp <- rast("exploration/colony summaries/crozet/thetao.nc")

# limit temporally
#temp <- temp[[time(temp) == as_date("2013-02-06")]]

# plot
test <- temp[[c(1, 25, 26, 27, 28, 29)]]
ggplot() + geom_spatraster(data = test) +
  geom_spatvector(data = sub_trax, size = 0.5) +
  facet_wrap(~lyr)

# plot for each date
dates <- sort(unique(time(temp)))
dates <- as_date("2013-02-05")
for(i in 1:length(dates)){
this_date <- dates[i]
slice <- temp[[time(temp) == this_date]]
slice <- slice[[25]]
p1 <- ggplot() + geom_spatraster(data = slice) +
  geom_spatvector(data = sub_trax, size = 0.5) +
  scale_fill_viridis_c() +
  ylim(-48.4, -46) +
  facet_wrap(~lyr) +
  ggtitle(this_date)
print(p1)
}

# plot ssh
ssh <- rast("E:/Satellite_Data/daily/ssh/ssh_2013.nc")

# limit temporally
ssh <- ssh[[time(ssh) == as_date("2013-02-06")]]
ssh <- crop(ssh, ext(trax) + c(1, 1, 1, 1))

# plot
ggplot() + geom_spatraster(data = ssh) +
  geom_spatvector(data = sub_trax, size = 0.5)

# plot uo and vo
uo <- rast("E:/Satellite_Data/daily/uo/uo_2013.nc")
vo <- rast("E:/Satellite_Data/daily/vo/vo_2013.nc")

# limit temporally
uo <- uo[[time(uo) == as_date("2013-02-06")]]
vo <- vo[[time(vo) == as_date("2013-02-06")]]
uo <- crop(uo, ext(trax) + c(1, 1, 1, 1))
vo <- crop(vo, ext(trax) + c(1, 1, 1, 1))

# calculate current
current <- sqrt(uo^2 + vo^2)

# plot
ggplot() + geom_spatraster(data = current) +
  geom_spatvector(data = sub_trax, size = 0.5)
