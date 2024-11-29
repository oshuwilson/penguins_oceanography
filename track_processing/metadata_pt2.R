#change timezone to UTC
zonediff <- readline("enter the time difference relative to UTC:")
tracks$datetime <- tracks$datetime - hours(zonediff)

#cleanup
rm(i, onefile, zonediff)


# 2. Read in Existing Metadata and keep relevant columns
#read in existing metadata
ext_meta <- read.csv(paste0("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Data/", species_code, "/", study_code, "/metadata.csv"))
