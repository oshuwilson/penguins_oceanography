
trip_split <- function(tracks, meta, buff.dist){

#get deployment sites
sites <- meta %>% 
  dplyr::select(individual_id, deployment_decimal_latitude, deployment_decimal_longitude)

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
buff <- buffer(colony, buff.dist)

#create plot for visualising buffer
p1 <- ggplot() + geom_spatvector(data = tracks, size = 0.5) +
  geom_spatvector(data = starts, col = "red", size = 2) + 
  geom_spatvector(data = buff, alpha = 0, lwd = 1, col = "red") +
  theme_bw()

#calculate distance to colony for each point
tracks <- tracks %>% 
  mutate(dist_to_colony = as.numeric(distance(tracks, colony)))

#convert back to dataframe
tracks <- as.data.frame(tracks, geom = "XY")

#identify points where distance to colony is less than 10000m
home_points <- tracks %>% 
  mutate(home = ifelse(dist_to_colony < buff.dist, 1, 0)) %>%
  filter(home == 1) %>%
  dplyr::select(-home)

#group home points into trip segments based on time differences of over 4 hours
home_points <- home_points %>% 
  arrange(individual_id, date) %>%
  group_by(individual_id) %>%
  mutate(lag_date = lag(date)) %>%
  mutate(timediff = difftime(date, lag_date, units = "hours")) %>%
  mutate(trip = ifelse(timediff > 4, 1, 0)) %>% 
  mutate(trip = ifelse(is.na(trip), 0, trip)) %>%
  mutate(trip = cumsum(trip)) %>%
  dplyr::select(-timediff, -lag_date)

#identify start/end points of trips (as the minimum distance to a colony within a timeframe)
terminals <- home_points %>% 
  group_by(individual_id, trip) %>%
  filter(dist_to_colony == min(dist_to_colony)) %>%
  ungroup() %>%
  dplyr::select(-trip) %>%
  mutate(trip_start = 1)

#append start_end points to tracks
tracks <- tracks %>%
  mutate(trip_start = 0) %>%
  anti_join(terminals, by = c("individual_id", "date", "dist_to_colony", "x", "y")) %>%
  bind_rows(terminals)

#group tracks into trips
tracks <- tracks %>%
  arrange(individual_id, date) %>%
  group_by(individual_id) %>%
  mutate(trip = cumsum(trip_start)) %>%
  mutate(trip = paste(individual_id, trip, sep = "_")) %>%
  ungroup()

#create time to start column
tracks <- tracks %>% 
  group_by(trip) %>%
  mutate(start_time = first(date)) %>%
  mutate(time_since_start = as.numeric(difftime(date, start_time, units = "hours"))) %>%
  dplyr::select(-start_time, -trip_start)

#remove trips that never leave boundary
short_trips <- tracks %>% 
  group_by(trip) %>%
  summarise(max_dist = max(dist_to_colony)) %>%
  filter(max_dist < buff.dist) %>%
  mutate(trip = as.factor(trip))
tracks <- tracks %>% 
  filter(!trip %in% levels(short_trips$trip)) %>%
  ungroup()

#remove trips with fewer than 10 points
small_trips <- tracks %>% 
  group_by(trip) %>%
  summarise(n = n()) %>%
  filter(n < 10) %>%
  mutate(trip = as.factor(trip))
tracks <- tracks %>%
  filter(!trip %in% levels(small_trips$trip)) %>%
  ungroup()

print(p1)
return(tracks)

}