#-------------------------------------------------------------------------------
# Remove obvious errors from HMM state assigment
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)

# define species, site, and stage
this.species <- "EMPE"
this.site <- "Taylor Glacier"
this.stage <- "post-breeding"

# read in state-assigned tracks
tracks <- readRDS(file = paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_unchecked.rds"))
nID <- length(unique(tracks$ID))
unique(as.character(tracks$ID))

# list trips for removal (using PDF plot)
# trips can be removed if 
# a) erroneous loops created by aniMotum become ARS
# b) trips are all one category (generally very short trips)
rm_trips <- c()

# remove trips with poor state assignments
tracks <- tracks %>%
  mutate(ID = as.character(ID)) %>%
  filter(!ID %in% rm_trips)

# check that IDs have been successfully removed
newn <- length(unique(tracks$ID))
if(nID - newn == length(rm_trips)){
  print("success")
} else {
  print("fail - please check all removal IDs are correct")
  rm(list=ls())
}

# export
saveRDS(tracks, file = paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_checked.rds"))

