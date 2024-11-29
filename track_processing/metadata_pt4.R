### deployment_site ###
#join from ext_meta if present
if(exists("deployment_site", where = ext_meta)){
  meta <- meta %>% 
    left_join(select(ext_meta, individual_id, deployment_site))
}

#if not present in ext_meta
if(!exists("deployment_site", where = meta)){
  meta$deployment_site <- readline("enter deployment site of all individuals if the same:")
}


### deployment date, time, and lat/lon ###
#join date/time if in ext_meta
if(exists("deployment_date", where = ext_meta)){
  meta <- meta %>% 
    left_join(select(ext_meta, individual_id, device_id, deployment_date)) %>%
    mutate(deployment_year = year(deployment_date),
           deployment_month = month(deployment_date),
           deployment_day = day(deployment_date)) %>% 
    select(-deployment_date)
}

if(exists("deployment_time", where = ext_meta)){
  meta <- meta %>% 
    left_join(select(ext_meta, individual_id, device_id, deployment_time)) %>%
    mutate(deployment_time = as_hms(deployment_time)) %>%
    mutate(deployment_time = as.character(deployment_time))
}


#join lat/lon if in ext_meta
if(exists("deployment_latitude", where = ext_meta)){
  meta <- meta %>% 
    left_join(select(ext_meta, individual_id, device_id, deployment_lat, deployment_lon))
}

#if date, time, and/or lat/lons are not in ext_meta, use first point for each device_ID, but exercise caution
if(!exists("deployment_date", where = ext_meta) | !exists("deployment_latitude", where = ext_meta)){
  tracks <- tracks %>%
    left_join(select(ext_meta, individual_id, device_id))
  
  deployments <- tracks %>% 
    group_by(individual_id, device_id) %>% 
    drop_na(all_of(c("lat", "lon"))) %>%
    summarise(start_time = first(datetime),
              start_lat = first(lat),
              start_lon = first(lon))
  
  #visualise to check
  deps <- vect(deployments, geom = c("start_lon", "start_lat"), crs = "EPSG:4326")
  deps <- project(deps, "EPSG:6932")
  e <- ext(deps) + c(10000, 10000, 10000, 10000)
  crop_coast <- crop(coast, e)
  plot(crop_coast)
  plot(deps, add=T)
  
  #extract deployment time and date
  deployments <- deployments %>% 
    mutate(deployment_year = year(start_time),
           deployment_month = month(start_time),
           deployment_day = day(start_time),
           deployment_time = as_hms(start_time)) %>%
    mutate(deployment_time = as.character(deployment_time))
  
  #format columns for metadata
  deployments <- deployments %>% 
    rename(deployment_decimal_latitude = start_lat,
           deployment_decimal_longitude = start_lon) %>%
    select(-start_time)
  
  #remove deployment date, time, or lat/lon columns if already in meta
  if(exists("deployment_day", where = meta)){
    deployments <- deployments %>% select(-deployment_year, -deployment_month, -deployment_day)
  }
  
  if(exists("deployment_time", where = meta)){
    deployments <- deployments %>% select(-deployment_time)
  }
  
  if(exists("deployment_decimal_latitude", where = meta)){
    deployments <- deployments %>% select(-deployment_decimal_latitude, -deployment_decimal_longitude)
  }
  
  #if happy, extract deployment information - only choose data missing from ext_meta (FIX THIS BIT)
  meta <- meta %>% 
    left_join(deployments,
              by = c("individual_id", "device_id"))
  
  #cleanup
  rm(deployments, deps, e, crop_coast)
}

### data contact and email ###
meta$data_contact <- readline("Enter Data Contact:")
meta$contact_email <- readline("Enter Contact Email:")


### split tracks into separate files and store filenames ###
#metadata filename column
meta$file_name <- files

### keepornot column for filtering ###
#NA for now - will be populated later
meta$keepornot <- NA


# 4. Join with existing metadata
#reorder columns
meta <- meta %>% select(dataset_identifier, file_name, individual_id, keepornot, scientific_name, 
                        common_name, abbreviated_name, sex, age_class, device_id, device_type, 
                        deployment_site, deployment_year, deployment_month, deployment_day,
                        deployment_time, deployment_decimal_longitude, deployment_decimal_latitude,
                        data_contact, contact_email)

#check whether there are duplicates before joining dataset
inds <- meta %>% 
  mutate(individual_id = as.character(individual_id)) %>%
  pull(individual_id)
devs <- meta %>% 
  mutate(device_id = as.character(device_id)) %>%
  pull(device_id)

dups <- all_meta %>% 
  filter(individual_id %in% inds & abbreviated_name == species_code |
           device_id %in% devs & abbreviated_name == species_code)

if(nrow(dups) > 0){
  print(dups)
  stop("Duplicate IDs found in metadata - check for temporal overlap in individuals or devices from dups dataset.")
}

#join together
meta <- meta %>% mutate(device_id = as.character(device_id))

all_meta <- all_meta %>% bind_rows(meta) %>% #remove duplicates
  arrange(abbreviated_name, dataset_identifier, individual_id) 

