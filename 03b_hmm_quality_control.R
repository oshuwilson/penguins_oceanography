#-----------------------------------
# Remove questionably assigned trips
#-----------------------------------

setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)

#clear environment
rm(list=ls())

# define species, site, and stage
this.species <- "CHPE"
this.site <- "South Shetland"
this.stage <- "post-moult"

# read in state-assigned tracks
tracks <- readRDS(file = paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_unchecked.rds"))
nID <- length(unique(tracks$ID))

# list trips for removal (using PDF plot)
# trips can be removed if 
# a) erroneous loops created by aniMotum become ARS
# b) ARS and transit behaviour is visually questionable
# c) trips are all one category (generally very short trips)
rm_trips <- c("165151_0", "165159_0", "165181_0", "165257_0", "CHPE188_0", "CHPE189_0", "CHPE199_0",
              "CHPE202_0", "CHPE203_0", "CHPE207_0", "CHPE216_0")

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
  rm(list=ls())
}

# export
saveRDS(tracks, file = paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_checked.rds"))

