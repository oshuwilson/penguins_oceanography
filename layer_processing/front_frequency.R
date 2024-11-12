#rolling front frequency algorithm

rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/")

{
  library(terra)
  library(lubridate)
  library(tidyterra)
}

#loop to process all years at once
for(z in 1993:2023) {
  try({
    
    #read in fronts for this year
    fronts <- readRDS(paste0("fronts/fronts_", z, ".RDS")) 
    
    #read fronts from previous year
    y <- z - 1
    try(old_fronts <- readRDS(paste0("fronts/fronts_", y, ".RDS"))) 
    
    #set day
    for(i in 1:nlyr(fronts)){
        
        #extract date of corresponding layer
        enddate <- time(fronts[[i]])
        
        #extract date 30 days before
        startdate <- enddate - days(30)
        
        #read in preceding thirty days - if all from same year
        if(i > 30){
          fr <- fronts[[time(fronts) >= startdate & time(fronts) <= enddate]]
        }
        
        #read in preceding thirty days - if from different years
        if(i <= 30){
          fr1 <- old_fronts[[time(old_fronts) >= startdate]]
          fr2 <- fronts[[time(fronts) <= enddate]]
          fr <- c(fr1, fr2)
          rm(fr1, fr2)
        }
        
        #front sums
        fsum <- sum(fr, na.rm=T)
        
        #sum number of days with data present for each cell
        nam <- c(NA, NA, 0,
                 -Inf, Inf, 1)
        nam <- matrix(nam, ncol=3, byrow=T)
        naclass <- classify(fr, nam)
        nasum <- sum(naclass)
        
        #front frequency
        freq <- fsum/nasum
        
        #append time info
        time(freq) <- enddate
        
        #if first slice in year create yearly front frequency SpatRaster
        if(!exists("frontfreq")){
          frontfreq <- freq
        }
        
        #if later slice, append to SpatRaster
        if(exists("frontfreq")){
          frontfreq <- c(frontfreq, freq)
        }
        
        #complete
        print(paste0(i, "/", nlyr(fronts)))
        
        #cleanup
        rm(fsum, freq, nasum, naclass, nam, fr)
    }
    
    #export front frequency file
    saveRDS(frontfreq, 
            file = paste0("front_freq/front_freq_", z, ".RDS"))
    
    #print when complete
    print(paste0(z, " completed"))
    
    #remove front frequency to avoid messing up future years
    rm(frontfreq)
    
  })
}