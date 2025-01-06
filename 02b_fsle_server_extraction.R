rm(list=ls())
setwd("/mnt/home/joshua/penguin_tracks")

library(raadtools)
library(terra)
library(lubridate)
library(dplyr)

#set species, site, and stage
this.species <- "KIPE"

#list all extraction files in this directory
files <- list.files(path = this.species, pattern = "*.RDS", full.names = T)

#for each file
for(file in files){
  
  #read in tracks
  tracks <- readRDS(file)
  
  #get list of x, y and time information from tracks
  a <- data.frame(x = tracks$x, y = tracks$y, time = tracks$date)
  
  #change negative lon values to 0-360 lon
  a <- a %>%
    mutate(x = ifelse(x<0, 360 + x, x))
  
  #extract fsle
  fsle <- raadtools::extract(readfsle, a)
  tracks$fsle <- fsle
  
  #list name of file without extension and folder
  clean <- basename(file)
  clean <- tools::file_path_sans_ext(clean)
  
  #export
  saveRDS(tracks, paste0("fsle/", this.species, "/", clean, "_with_fsle.RDS"))
  
  #print to signify completion
  print(clean)
  
}
