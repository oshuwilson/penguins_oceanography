# 3. Create dataframe for individual study

### group by individual_id and device_id, and append study_code and species_code ###
meta <- tracks %>% group_by(individual_id, device_id) %>%
  summarise(dataset_identifier = study_code, abbreviated_name = species_code) %>%
  ungroup()

### scientific_name, common_name ### 
#import list of names
names <- read.csv("~/OneDrive - University of Southampton/Documents/RAATD 2.0/Metadata/species_codes.csv")

#filter to this species
names <- filter(names, abbreviated_name == species_code)

#assign common and scientific name from this
meta$scientific_name <- names$scientific_name
meta$common_name <- names$common_name

#cleanup
rm(names)


### sex, age_class ###
#sex - if in ext_meta
if(exists("sex", where = ext_meta)){
  meta <- meta %>% 
    left_join(select(ext_meta, individual_id, sex))
}

#sex - if not in ext_meta (if different for different individuals do this manually)
if(!exists("sex", where = ext_meta)){
  meta$sex <- readline("enter sex of all individuals if the same (male, female, or unknown):")
  
  if(first(meta$sex) == "NA"){
    meta$sex <- NA
  }
}

#check that sex is either male, female, or NA
poss_sexes <- c("male", "female", "unknown")
invalid <- meta %>%
  filter(!sex %in% poss_sexes)
if(nrow(invalid) > 0){
  print(invalid)
  stop("sex must be male, female or unknown. Please rename others.")
}

#cleanup
rm(invalid, poss_sexes)

#age_class - if in ext_meta
if(exists("age_class", where = ext_meta)){
  meta <- meta %>% 
    left_join(select(ext_meta, individual_id, age_class))
}

#age_class - if not in ext_meta (define as Adult or Juvenile if known and NA if not)
if(!exists("age_class", where = ext_meta)){
  meta$age_class <- readline("enter age class of all individuals if the same (adult, juvenile, or unknown):")
  
  if(first(meta$age_class) == "NA"){
    meta$age_class <- NA
  }
}

#check that age_class is either Adult, Juvenile, or NA
poss_ages <- c("adult", "juvenile", "unknown")
invalid <- meta %>%
  filter(!age_class %in% poss_ages)
if(nrow(invalid) > 0){
  print(invalid)
  stop("age_class must be adult, juvenile, or unknown. Please rename others.")
}

#cleanup
rm(invalid, poss_ages)


### device_type ###

#device_type - if in ext_meta
if(exists("device_type", where = ext_meta)){
  meta <- meta %>% 
    left_join(select(ext_meta, individual_id, device_type))
}

#device_type - if not in ext_meta (define as GPS, PTT, or GLS)
if(!exists("device_type", where = ext_meta)){
  meta$device_type <- readline("enter device type of all individuals if the same (GPS, PTT, or GLS):")
}

#check that device_type has been saved as GPS, PTT, or GLS
poss_devices <- c("GPS", "PTT", "GLS")
invalid <- meta %>%
  filter(!device_type %in% poss_devices)
if(nrow(invalid) > 0){
  print(invalid)
  stop("device_type must be GPS, PTT, or GLS. Please rename others.")
}

#save copy to prevent full restart
meta3_banked <- meta

#cleanup
rm(invalid, poss_devices)
