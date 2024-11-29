#formatting RAATD metadata - run from source (Ctrl+Shift+S) NEED DIFFERENT META NAME FOR EACH STEP
{
  library(tidyverse)
  library(tidyterra)
  library(hms)
}

# 1. Read in and format track dataframe

#set species, study and colony/region code
species_code <- readline("enter the 4-letter RAATD species code:") 
study_code <- readline("enter the dataset identifier (the folder name):") 

#check that study hasn't already been processed before beginning
all_meta <- readRDS("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/RAATD_2_Metadata.RDS")
if(study_code %in% all_meta$dataset_identifier){
  reps <- all_meta %>% filter(dataset_identifier == study_code, abbreviated_name == species_code)
  if(nrow(reps) > 0){
  stop("The metadata for this study has already been processed.")
  }
  rm(reps)
}

#set working directory to read files from
setwd(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code, "/Data"))

#list all available files
files <- list.files(pattern="*.csv") #change pattern depending on filetype

#are all tracks in one file?
onefile <- if_else(length(files) == 1, T, F)

#read tracks - if all IDs in one file
if(onefile == TRUE){
  tracks <- read.csv(files)
}

#read tracks - if IDs are in separate files
if(onefile == FALSE){
  first_num <- readline("enter the number of characters from the start of filenames where individual_id begins:")
  last_num <- readline("enter the number of characters from the end of filenames where individual_id ends:")
  tracks <- NULL
  for(i in files){
    fixes <- read.csv(i)
    fixes$individual_id <- substr(i, first_num, last_num) 
    tracks <- bind_rows(tracks, fixes)
    rm(fixes)
  }
}
