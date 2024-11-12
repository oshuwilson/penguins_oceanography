#daily distance to ice edge
rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/")

{
  library(terra)
  library(sf)
  library(tidyterra)
  library(lubridate)
  library(ncdf4)
}

#ice sheets, tongues, and rumples
perm_ice <- readRDS("coast/ice_features.RDS")

#loop to process all years
for(z in 1993:2023) {
    
    rm(dist2ice)
    
    #read in yearly sic
    sic <- rast(paste0("sic/sic_", z, ".nc"))
    
    #run loop over each day
    for(i in 1:nlyr(sic)){
    
    try({
      #remove variables from previous loop runs 
      rm(list=setdiff(ls(), c("perm_ice", "z", "i", "sic", "dist2ice")))
      
      #extract SIC for that day
      sic_daily <- sic[[i]]
      
      #crop to area of interest 
      e <- ext(-180, 180, -80, -40)
      sic_daily <- crop(sic_daily, e)
      
      #classify SIC values as <0.15 or >0.15
      m <- c(-Inf, 0.15, NA,
             0.15, Inf, 0.15)
      m <- matrix(m, ncol=3, byrow = T)
      classes <- classify(sic_daily, m)
      
      #convert to polygons
      poly <- as.polygons(classes)
      
      #project to Antarctic view
      poly <- project(poly, crs(perm_ice))
      
      #merge into one object and fill gaps within sea ice edge
      outline <- rbind(perm_ice, poly)
      outline <- aggregate(outline)
      outline <- fillHoles(outline)
      
      #rasterize to stereographic raster
      depth <- rast("depth/depth_stereographic.RDS")
      ice_edge <- rasterize(outline, depth)
      
      #calculate distance to ice edge
      dist <- distance(ice_edge)
      
      #invert and calculate distance to ice edge within sea ice extent (for polynyas)
      m2 <- c(NA, NA, 1,
              0, Inf, 0)
      m2 <- matrix(m2, ncol=3, byrow=T)
      open_water <- classify(ice_edge, m2)
      
      inverted <- distance(open_water, target=0)
      inverted <- inverted * -1
      
      #combine two rasters
      final_dist <- dist + inverted
      
      #round to nearest km
      final_dist <- final_dist/1000
      final_dist <- round(final_dist)
      
      #project to original SIC crs
      dist2iceedge <- project(final_dist, crs(sic_daily))
      
      #recrop to area of interest
      e <- ext(-180, 180, -80, -40)
      dist2iceedge <- crop(dist2iceedge, e)
      
      #resample to GLORYS resolution
      dist2iceedge <- resample(dist2iceedge, sic_daily, method = "bilinear")
      
      #assign time to layer
      time(dist2iceedge) <- time(sic_daily)
      
      #combine layer into yearly raster
      if(i == 1){
        dist2ice <- dist2iceedge
      }
      
      if(i != 1){
        dist2ice <- c(dist2ice, dist2iceedge)
      }
      
      #print completion
      print(paste0(i, "/", nlyr(sic)))
      
      #remove temporary files
      tmpFiles(remove = TRUE)
      
      })
    }
    
    #export yearly raster
    writeCDF(dist2ice,
             filename = paste0("dist2ice/dist2ice_", z, ".nc"),
             varname = "dist2ice",
             longname = "distance to 15% sea ice edge",
             unit = "km",
             overwrite = TRUE,
             prec = "integer") #export as INT2S to minimise memory use
    
    #cleanup to avoid impacting future years
    rm(dist2ice, sic)
    
    #print completion
    print(paste0(z, " completed"))
  
}
