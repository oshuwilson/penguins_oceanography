#-----------------------------------
# remove questionably assigned trips
#-----------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)

# define species, site, and stage
this.species <- "KIPE"
this.site <- "Marion"
this.stage <- "chick-rearing"

# read in state-assigned tracks
tracks <- readRDS(file = paste0("output/hmm/hmm_tracks/", this.species, "_", this.site, "_", this.stage, "_tracks_unchecked.rds"))

# list trips for removal (using PDF plot)
# trips can be removed if 
# a) erroneous loops created by aniMotum become ARS
# b) ARS and transit behaviour is visually questionable
# c) trips are all one category (generally very short trips)
rm_trips <- c("57335_4", "57335_1", "119308_5", "119305_1",
              "57332_1", "57331_0", "119308_1",
              "119312_0", "119309_3", "119305_2")

# remove trips with poor state assignments
tracks <- tracks %>%
  filter(!ID %in% rm_trips)

# export
saveRDS(tracks, file = paste0("output/hmm/hmm_tracks/", this.species, "_", this.site, "_", this.stage, "_tracks_checked.rds"))

