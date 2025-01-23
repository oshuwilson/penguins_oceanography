#-----------------------------------
# Remove questionably assigned trips
#-----------------------------------

setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)

#clear environment
rm(list=ls())

# define species, site, and stage
this.species <- "MAPE"
this.site <- "South Georgia"
this.stage <- "post-breeding"

# read in state-assigned tracks
tracks <- readRDS(file = paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_unchecked.rds"))
nID <- length(unique(tracks$ID))

# list trips for removal (using PDF plot)
# trips can be removed if 
# a) erroneous loops created by aniMotum become ARS
# b) ARS and transit behaviour is visually questionable
# c) trips are all one category (generally very short trips)
rm_trips <- c("MAPE-dtsetBirdLife748-713_5-RAATD_1", "MAPE-dtsetBirdLife748-713_5-RAATD_2", "MAPE-dtsetBirdLife748-713_5-RAATD_3",
              "MAPE-dtsetBirdLife748-713_5-RAATD_4", "MAPE-dtsetBirdLife748-713_7-RAATD_0", "MAPE-dtsetBirdLife751-P1-RAATD_2",
              "MAPE-dtsetBirdLife751-P2-RAATD_1", "MAPE-dtsetBirdLife751-P6-RAATD_1", "MAPE-dtsetBirdLife751-X32-RAATD_2",
              "MAPE-dtsetBirdLife751-depid483-RAATD_1")

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

