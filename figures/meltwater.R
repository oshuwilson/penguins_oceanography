#--------------------
# Meltwater Influence
#--------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(marmap)
}

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/EMPE/Rothschild Island post-breeding tracks checked.rds")
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")


# find individuals of interest
ed <- trax %>%
  filter(eddies <= -1 | eddies >= 1) %>%
  filter(state == "ARS")
plot(ed)
plot(trax, add = T, col = "red", pch = ".")

ggplot(trax, aes(col = individual_id)) + geom_spatvector()

# individuals of interest:
# - 98534
# - 131004
# - 150624
# - 131432
# - 98532

# filter to these individuals
trax <- trax %>%
  filter(individual_id %in% c(98534, 131004, 150624, 131432, 98532))

# read in SIC for 2015 and 2016
sic1 <- rast("E:/Satellite_Data/daily/sic/sic_2015.nc")
sic2 <- rast("E:/Satellite_Data/daily/sic/sic_2016.nc")

# limit to dates of tracks
dates <- unique(as_date(trax$date))
sic1 <- sic1[[time(sic1) %in% dates]]
sic2 <- sic2[[time(sic2) %in% dates]]

# crop to track extent
e <- ext(trax) + c(1, 1, 1, 1)
sic1 <- crop(sic1, e)
sic2 <- crop(sic2, e)

# join sic rasters together
sic <- c(sic1, sic2)

# for each individual, identify eddy foraging hotspots and plot
inds <- unique(trax$individual_id)
ind <- inds[5]

# limit eddy foraging points to this individual
this_ed <- ed %>%
  filter(individual_id == ind)

# for each date in foraging events, plot SIC and tracks
ed_dates <- unique(as_date(this_ed$date))
for(i in 1:length(ed_dates)){
  this_date <- ed_dates[i]
  this_sic <- sic[[time(sic) == this_date]]
  this_ed_trax <- this_ed %>%
    filter(as_date(date) == this_date)
  p1 <- ggplot() + geom_spatraster(data = this_sic) +
    geom_spatvector(data = this_ed_trax, size = 2, col = "darkred") + 
    ggtitle(paste0(ind, " ", this_date))
  print(p1)
}


# dates of interest

# IND 131004
# - 3rd to 4th December 2015

# IND 98532
# - 13th to 16th April 2016

# IND 98534
# - 17th to 22nd Jan 2016


# 1. Individual 131004 3rd to 4th December 2015

# filter trax to this individual
this_trax <- trax %>%
  filter(individual_id == 131004)

# filter sic to these dates
this_sic <- sic[[time(sic) %in% as_date(c("2015-12-03", "2015-12-04"))]]

# filter trax to between these dates (+/- 2 days)
this_trax <- this_trax %>%
  filter(as_date(date) >= as_date("2015-12-01") & as_date(date) <= as_date("2015-12-06"))

# crop sic to extent of trax (plus 1 degree)
e2 <- ext(this_trax) + c(1, 1, 1, 1)
this_sic <- crop(this_sic, e2)

# plot
ggplot() + geom_spatraster(data = this_sic) +
  geom_spatvector(data = this_trax, size = 2, aes(col = state)) +
  facet_wrap(~lyr)


# 2. Individual 98532 13th to 16th April 2016

# filter trax to this individual
this_trax <- trax %>%
    filter(individual_id == 98532)

# filter sic to these dates
this_sic <- sic[[time(sic) %in% as_date(c("2016-04-13", "2016-04-14", "2016-04-15", "2016-04-16"))]]

# filter trax to between these dates (+/- 2 days)
this_trax <- this_trax %>%
    filter(as_date(date) >= as_date("2016-04-11") & as_date(date) <= as_date("2016-04-18"))

# crop sic to extent of trax (plus 1 degree)
e2 <- ext(this_trax) + c(1, 1, 1, 1)
this_sic <- crop(this_sic, e2)

# plot
ggplot() + geom_spatraster(data = this_sic) +
    geom_spatvector(data = this_trax, size = 2, aes(col = state)) +
    facet_wrap(~lyr)


# 3. Individual 98534 17th to 22nd Jan 2016

# filter trax to this individual
this_trax <- trax %>%
    filter(individual_id == 98534)

# filter sic to these dates
this_sic <- sic[[time(sic) %in% as_date(c("2016-01-17", "2016-01-18", "2016-01-19", "2016-01-20", "2016-01-21", "2016-01-22"))]]

# filter trax to between these dates (+/- 2 days)
this_trax <- this_trax %>%
    filter(as_date(date) >= as_date("2016-01-15") & as_date(date) <= as_date("2016-01-24"))

# crop sic to extent of trax (plus 1 degree)
e2 <- ext(this_trax) + c(1, 1, 1, 1)
this_sic <- crop(this_sic, e2)

# plot
ggplot() + geom_spatraster(data = this_sic) +
    geom_spatvector(data = this_trax, size = 2, aes(col = state)) +
    facet_wrap(~lyr)


# Plot salinity for #1
# filter trax to this individual
this_trax <- trax %>%
    filter(individual_id == 131004)

# read in salinity
sal <- rast("E:/Satellite_Data/daily/sal/sal_2015.nc")
sal <- sal[[time(sal) %in% as_date(c("2015-12-03", "2015-12-04"))]]

# crop salinity to track extent
sal <- crop(sal, ext(this_trax))

# plot
ggplot() + geom_spatraster(data = sal) +
  geom_spatvector(data = this_trax, size = 2, aes(col = state)) +
  scale_fill_viridis_c() +
  ylim(-67, -65) +
  xlim(-81, -78) +
  facet_wrap(~lyr)


# Plot ssh for #1

# read in ssh
ssh <- rast("E:/Satellite_Data/daily/ssh/ssh_2015.nc")
ssh <- ssh[[time(ssh) %in% as_date(c("2015-12-03", "2015-12-04"))]]

# crop ssh to track extent
ssh <- crop(ssh, ext(this_trax))

# plot
ggplot() + geom_spatraster(data = ssh) +
  geom_spatvector(data = this_trax, size = 2, aes(col = state)) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr)

# plot ssh for #2

# filter trax to this individual
this_trax <- trax %>%
    filter(individual_id == 98532)

# filter tracks to dates of interest
this_trax <- this_trax %>%
    filter(as_date(date) >= as_date("2016-04-11") & as_date(date) <= as_date("2016-04-18"))

# read in ssh
ssh <- rast("E:/Satellite_Data/daily/ssh/ssh_2016.nc")
ssh <- ssh[[time(ssh) %in% as_date(c("2016-04-13", "2016-04-14", "2016-04-15", "2016-04-16"))]]

# crop ssh to track extent
ssh <- crop(ssh, ext(this_trax))

# plot
ggplot() + geom_spatraster(data = ssh) +
  geom_spatvector(data = this_trax, size = 2, aes(col = state)) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr)


# plot ssh for #3

# filter trax to this individual
this_trax <- trax %>%
    filter(individual_id == 98534)

# filter tracks to dates of interest
this_trax <- this_trax %>%
    filter(as_date(date) >= as_date("2016-01-15") & as_date(date) <= as_date("2016-01-24"))

# read in ssh
ssh <- rast("E:/Satellite_Data/daily/ssh/ssh_2016.nc")
ssh <- ssh[[time(ssh) %in% as_date(c("2016-01-17", "2016-01-18", "2016-01-19", "2016-01-20", "2016-01-21", "2016-01-22"))]]

# crop ssh to track extent
ssh <- crop(ssh, ext(this_trax))

# plot
ggplot() + geom_spatraster(data = ssh) +
  geom_spatvector(data = this_trax, size = 2, aes(col = state)) +
  scale_fill_viridis_c() +
  facet_wrap(~lyr)
