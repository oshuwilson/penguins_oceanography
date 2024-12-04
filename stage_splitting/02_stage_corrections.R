#--------------------------------------------------------------
# Correct breeding stages based on visual checks and literature
#--------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
}


#----------------
# Adelie Penguins
#----------------

this.species <- "ADPE"

#read in stages and tracks with stage trips
stages <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_stages.RDS"))
tracks <- readRDS(paste0("code/stage_splitting/stage_dates/", this.species, "_tracks_with_stage_trips.RDS"))

#append study codes from metadata
meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
meta <- meta %>% filter(abbreviated_name == this.species)


# 1. IAU_KG are all chick-rearing - see Machado-Gaye et al. (2024): doi.org/10.1007/s00227-024-04390-w

#correct iau_kgs to chick-rearing
iau_kg <- meta %>% filter(dataset_identifier == "IAU_KG")

#rename all stages to chick to chick-rearing
iau_kg_stages <- stages %>% 
  filter(individual_id %in% iau_kg$individual_id) %>%
  mutate(stage = "chick-rearing")

#if multiple stages present, choose earliest start date and latest end date
iau_kg_stages <- iau_kg_stages %>%
  group_by(individual_id, device_id) %>%
  summarise(start = min(start), end = max(end), stage = "chick-rearing")

#export
#saveRDS(iau_kg_stages, paste0("code/stage_splitting/stage_dates/", this.species, "_stages.RDS"))
