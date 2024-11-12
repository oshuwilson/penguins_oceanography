#daily front detection algorithm
rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/")

{
  library(terra)
  library(grec)
  library(foreach)
  library(doParallel)
}

#setup parallel programming
registerDoParallel(cores = 31)

#loop to process all years at once
foreach(z = 1993:2023) %dopar% {
  try({
    
    #read in ssh data
    ssh_year <- rast(paste0("ssh/ssh_", z, ".nc"))
    
    #run loop over each day
    for(i in 1:nlyr(ssh_year)){ 
      ssh <- ssh_year[[i]]
      
      #crop to area of interest - SMALL FOR NOW
      e <- ext(-180, 180, -80, -30)
      ssh <- crop(ssh, e)
      
      #calculate ssh gradients using Belkin O'Reilly contextual median filter
      hgrad <- getGradients(ssh, method = "BelkinOReilly2009")
      
      #identify fronts using 0.25m per km threshold recommended by Sokolov and Rintoul (2007)
      m <- c(-Inf, 0.25, 0,
             0.25, Inf, 1)
      m <- matrix(m, ncol = 3, byrow = TRUE)
      class <- classify(hgrad, m)
      plot(class)
      
      #append time info
      time(class) <- time(ssh)
      
      #if first slice in year create yearly front SpatRaster
      if(i == 1){
        fronts <- class
      }
      
      #if later slice, append to SpatRaster
      if(i > 1){
        fronts <- c(fronts, class)
      }
    }
    
    #export front file
    saveRDS(fronts, 
            file = paste0("fronts/fronts_", z, ".RDS"))
    
    #print when complete
    print(paste0(z, " completed"))
    
  })
}