# get test statistics from GAMMs

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(gamm4)
}

# read in species-site-stage metadata
srs <- read.csv("data/tracks/species_site_stage_v2.csv")

# remove sites where eddies were not encountered (no model created)
srs <- srs %>%
  filter(!site %in% c("Tierra del Fuego", "Cape Bird", "Mamejima Island", "Cape Colbeck")) %>%
  filter(site != "Edmonson Point" | stage == "pre-moult") %>%
  filter(site != "Point Thomas, King George Island" | stage != "chick-rearing") %>%
  filter(site != "Stranger Point, King George Island" | stage != "chick-rearing") %>%
  filter(site != "Signy" | stage != "incubation") %>%
  filter(site != "Pointe Geologie" | stage != "post-breeding") %>%
  filter(site != "Rothschild Island" | stage != "chick-rearing")

for(i in 38:nrow(srs)){
  
  # define species, site, and stage
  this.species <- srs$species[i]
  this.site <- srs$site[i]
  this.stage <- srs$stage[i]
  area <- srs$island[i]
  auger <- srs$auger[i]
  
  # read in gamm4 model for this case study
  if(auger == "yes"){
    m1 <- readRDS(paste0("output/models/", this.species, "/", this.species, " ", this.site, " ", this.stage, " auger gamm.rds"))
  } else {
    m1 <- readRDS(paste0("output/models/", this.species, "/", this.site, " ", this.stage, " gamm.rds"))
  }
  
  #-------------------------------------------------------
  # Get Key Criteria
  #-------------------------------------------------------
  # need the smooth term info (Chi-squared + p-value)
  # also the adjusted R-squared of the model
  
  # create the summary object
  summ <- summary(m1$gam)
  
  # get the smooth term values
  smooths <- summ$s.table %>%
    as.data.frame()
  
  # create variable column
  smooths$variable <- rownames(smooths)
  
  # change variable names to remove s()
  smooths$variable <- gsub("s\\(|\\)", "", smooths$variable)
  
  # get the adjusted R-squared for the model
  rsq <- summ$r.sq
  
  # add to table
  smooths$modelrsq <- rsq
  
  # add species site and stage to the table
  smooths$species <- this.species
  smooths$site <- this.site
  smooths$stage <- this.stage
  
  # join to all other case studies
  if(i == 1){
    all_smooths <- smooths
  } else {
    all_smooths <- rbind(all_smooths, smooths)
  }
  
  # print success
  print(paste0(this.species, " ", this.site, " ", this.stage, " complete"))
  
}


# export
saveRDS(all_smooths, "output/models/model_stats.rds")
write_csv(all_smooths, "output/models/model_stats.csv")
