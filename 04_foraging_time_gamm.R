#----------------------------------
# Assess eddy usage by feeding time 
#----------------------------------

rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/Chapter_02/")

{
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(gamm4)
  library(gratia)
  library(suncalc)
}

# read in species site stage info to loop over
srs <- read.csv("data/tracks/species_site_stage_v2.csv")

# sites of interest
# srs <- srs %>%
#   mutate(case_study = paste(species, site, stage, sep = " ")) %>%
#   filter(case_study %in% c("ADPE Windmill Islands chick-rearing",
#                            "ADPE Bechervaise Island chick-rearing",
#                            "ADPE Admiralty Bay, South Shetland chick-rearing",
#                            "ADPE Admiralty Bay, South Shetland incubation",
#                            "EMPE Taylor Glacier fledglings",
#                            "EMPE Auster Rookery fledglings",
#                            "ADPE Cape Bird chick-rearing")) %>%
#   select(-case_study)

# isolate colony and breeding stage
for(i in 1:nrow(srs)){
  this.species <- srs$species[i]
  this.site <- srs$site[i]
  this.stage <- srs$stage[i]
  area <- srs$island[i]
  
  # load in hmm checked tracks
  tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))
  
  # if number of distinct individuals is less than 3, skip
  if(n_distinct(tracks$individual_id) < 3){
    next
  }
  
  # 1. Process Data
  
  #read in original tracks to get lat/lons and error info
  original <- readRDS(paste0("output/tracks/", this.species, "/", area, " ", this.stage, " tracks.RDS"))
  
  #append latitudes, longitudes, and errors to state tracks
  tracks <- tracks %>% 
    left_join(select(original, individual_id, date, lon, lat, 
                     longitude_se, latitude_se, 
                     lon_se_km, lat_se_km))
  
  #remove tracks with large error
  tracks <- tracks %>%
    filter((latitude_se < 0.2 & longitude_se < 0.5 | #greater allowance for longitude as this can be compressed at poles
             (lon_se_km < 20000 & lat_se_km < 20000)))
  
  #create column in date format for suncalc
  tracks <- tracks %>%
    rename(datetime = date) %>%
    mutate(date = as_date(datetime))
  
  #get dawn times
  tracks$dawn <- getSunlightTimes(data = tracks,
                                  keep = c("dawn"), tz = "UTC") %>%
    pull(dawn)
  
  #get dusk times
  tracks$dusk <- getSunlightTimes(data = tracks,
                                  keep = c("dusk"), tz = "UTC") %>%
    pull(dusk)
  
  #only keep points between sunrise and sunset
  tracks <- tracks %>%
    filter(datetime >= dawn & datetime <= dusk |
             is.na(dawn) & is.na(dusk))
  
  # filter to ARS only
  ars <- tracks %>% filter(state == "ARS")
  
  # resample non-eddies to 0
  ars <- ars %>%
    mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))
  
  # read in background samples
  back <- readRDS(paste0("output/background/", this.species, "/", this.site, " ", this.stage, " background.rds"))
  
  # resample non-eddies to 0
  back <- back %>%
    mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))
  
  # create binary presence/absence cols
  ars$pa <- 1
  back$pa <- 0
  
  # select key variables
  ars <- ars %>% 
    select(individual_id, ed2, depth, curr, sic, pa)
  back <- back %>%
    select(individual_id, ed2, depth, curr, sic, pa)
  
  
  # 2. GAMM
  
  # join datasets together
  data <- bind_rows(ars, back)
  
  # remove points with SIC > 0.1
  data <- data %>%
    filter(sic <= 0.1 | is.na(sic)) %>%
    select(-sic)
  
  # if eddies vary in only one or no individuals, export smooth file as 0 and skip
  smallvar <- ars %>%
    select(ed2, individual_id) %>%
    group_by(individual_id) %>%
    summarise_all(n_distinct)
  ee <- smallvar %>% 
    filter(ed2 == 1) %>%
    nrow()
  if(ee >= nrow(smallvar) - 1){
    sm <- data.frame(ed2 = 0, .estimate = 0, .lower_ci = 0, .upper_ci = 0)
    write.csv(sm, paste0("output/gamms/smooths/", this.species, "/", this.site, " ", this.stage, " ed2.csv"))
    next
  }
  
  # run GAMM 
  m1 <- gamm4(pa ~ s(ed2, bs = "ts") + s(depth, bs = "ts") + s(curr, bs = "ts"),
              random = ~(1|individual_id), family = binomial, data = data)
  
  
  # 3. Odds Ratios
  
  # get smooths
  sm <- smooth_estimates(m1$gam, n = 1000) %>%
    add_confint()
  
  # apply exponential to smooths for odds ratios
  sm <- sm %>%
    mutate(.estimate = exp(.estimate),
           .lower_ci = exp(.lower_ci),
           .upper_ci = exp(.upper_ci))
  
  # plot
  p1 <- ggplot(sm) + 
    geom_line(aes(x = ed2, y = .estimate)) + 
    geom_ribbon(aes(x = ed2, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) + 
    theme_minimal() + 
    geom_hline(yintercept = 1, linetype = "dashed", col = "darkred") +
    labs(x = "Relative Eddy Distance", y = "Odds Ratio") +
    ylim(0, 4)
  
  
  # 4. Export
  
  # export model
  saveRDS(m1, paste0("output/gamms/models/", this.species, "/", this.site, " ", this.stage, " gamm.rds"))
  
  # export smooths
  saveRDS(sm, paste0("output/gamms/smooths/", this.species, "/", this.site, " ", this.stage, " smooths.rds"))
  
  # export gamm plot
  ggsave(paste0("output/gamms/plots/", this.species, "/", this.site, " ", this.stage, " gamm.png"), 
         p1, width = 10, height = 6, create.dir = T)
  
  # print completion
  print(paste0(this.species, " ", this.site, " ", this.stage, " complete"))
  
}

