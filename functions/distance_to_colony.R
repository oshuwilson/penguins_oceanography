distance_to_colony <- function(steps, meta){
  
#convert steps to terra
steps <- vect(steps, geom = c("x2_", "y2_"), crs = "epsg:6932")

#get deployment sites
sites <- meta %>% 
  select(individual_id, deployment_decimal_latitude, deployment_decimal_longitude)

#plot start points
starts <- vect(sites, 
               geom = c("deployment_decimal_longitude", "deployment_decimal_latitude"), 
               crs = "epsg:4326")
starts <- project(starts, "epsg:6932")

#if content, choose the average start point
colony_lon <- median(sites$deployment_decimal_longitude)
colony_lat <- median(sites$deployment_decimal_latitude)
colony <- as.data.frame(cbind(colony_lon, colony_lat))

#create buffer around colony - experiment with distance to decide on trip distance
colony <- vect(colony, geom = c("colony_lon", "colony_lat"), crs = "epsg:4326")
colony <- project(colony, "epsg:6932")

#create plot for visualising buffer
p1 <- ggplot() + geom_spatvector(data = steps, size = 0.5) +
  geom_spatvector(data = starts, col = "red", size = 2) +
  theme_bw()

#calculate distance to colony for each point
steps <- steps %>% 
  mutate(dist_to_colony = as.numeric(distance(steps, colony)))

#convert back to dataframe
steps <- as.data.frame(steps, geom = "XY")

print(p1)
return(steps$dist_to_colony)

}