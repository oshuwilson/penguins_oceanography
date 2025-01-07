#-----------------------------------
# Remove questionably assigned trips
#-----------------------------------

setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)

#clear environment
rm(list=ls())

# define species, site, and stage
this.species <- "KIPE"
this.site <- "Tierra del Fuego"
this.stage <- "incubation"

# read in state-assigned tracks
tracks <- readRDS(file = paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_unchecked.rds"))
nID <- length(unique(tracks$ID))

# list trips for removal (using PDF plot)
# trips can be removed if 
# a) erroneous loops created by aniMotum become ARS
# b) ARS and transit behaviour is visually questionable
# c) trips are all one category (generally very short trips)
rm_trips <- c("14A0001_14_2", "14A0487_14_1", "14A0490_14_2", "14A0491_16_4", 
              "153442_16._6", "153442_16._7", "153444_16._10")

# remove trips with poor state assignments
tracks <- tracks %>%
  mutate(ID = as.character(ID)) %>%
  filter(!ID %in% rm_trips)

# check that IDs have been successfully removed
newn <- length(unique(tracks$ID))
if(nID - newn == length(rm_trips)){
  print("success")
} else {
  print("fail")
}

# export
saveRDS(tracks, file = paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_checked.rds"))

