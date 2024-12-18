#----------------------------
# Create yearly polynya files
#----------------------------

#EPSG3412

rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/")

library(terra)
library(dplyr)
library(lubridate)
library(tidyterra)
library(ncdf4)

#list all .nc files 
files <- list.files(path = "polynyas/", pattern = ".nc", full.names = T)

#remove lonlat from list
files <- files[-4078]

#load in longitudes and latitudes
lonlat <- nc_open("polynyas/LonLat.nc")
lon <- ncvar_get(lonlat, "Lon")
lat <- ncvar_get(lonlat, "Lat")

#loop over every year
for(i in 2020:2022){
  
  #define year
  this.year <- i
  
  #find all files from this year
  yearly_files <- files[substr(files, 31, 34) == this.year]
  
  #for each day
  for(file in yearly_files){
  
  try({
      
    poly <- nc_open(file)
    
    #get values
    map <- ncvar_get(poly, "Map")
    
    #convert lon, lat, and values to cols to form dataframe
    lon <- as.vector(lon)
    lat <- as.vector(lat)
    map <- as.vector(map)
    
    #create matrix for rasterisation
    df <- data.frame(lon, lat, map)
    
    #replace values above 60 degrees south with NA
    df <- df %>% mutate(map = ifelse(lat < -60.0, map, NA)) 
    
    #reclassify map values to binary code
    df <- df %>% mutate(map = ifelse(map > 0, 1, 0))
    
    #create raster of polynyas in native projection
    pol <- rast(ncols = 1328, nrows = 1264,
                xmin = -3325000, xmax = 3325000, #epsg 3412 x-limits
                ymin = -3325000, ymax = 3325000, #espg 3412 y-limits
                crs = "epsg:3412")
    values(pol) <- df$map
    e <- ext(pol)
    
    #remove empty space and redefine co-ordinate limits (to match projection 3412)
    pol <- trim(pol, padding = 0)
    ext(pol) <- e 
    
    #rotate 90 degrees to remove longitudinal offset
    pol <- t(flip(pol, "horizontal"))
    
    #assign date from file name
    date <- substr(file, 31, 38)
    date <- as_date(date, format = "%Y%m%d")
    time(pol) <- date
    
    #stack with other rasters from the same year
    if(file == yearly_files[1]){
      pol_stack <- pol
    } else {
      pol_stack <- c(pol_stack, pol)
    }
    
    #print completion
    print(paste0(date, " complete"))
    
    })
    
  }
  
  #export
  writeCDF(pol_stack, file = paste0("polynyas//processed//polynyas_", this.year, ".nc"), overwrite = T)
  
  #print completion
  print(paste0(this.year, " complete"))
  
}
