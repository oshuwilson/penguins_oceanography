#----------------------------------
# Assess eddy usage by feeding time 
#----------------------------------

# changes to apply to SIC GAMMs as well
# 1. reduce error threshold (done here)
# 2. standardise covariates (done here)
# 3. account for spatiotemporal autocorrelation (separate preprocessing script, gamm correlation check here)
# 4. collect GAMM prediction metrics (e.g. AUC, CBI) from models excluding some individuals (done here) 
# 5. fit a GAMM without eddies and collect AIC/BIC statistics for model selection (done here)

rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/Chapter_02/")

{
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(terra)
  library(gamm4)
  library(gratia)
  library(suncalc)
  library(tidysdm)
}

# read in species site stage info to loop over
srs <- read.csv("data/tracks/species_site_stage_v2.csv")
srs <- srs %>% filter(species %in% c("KIPE", "MAPE"))

# sites of interest

# isolate colony and breeding stage
for(i in 18){
  this.species <- srs$species[i]
  this.site <- srs$site[i]
  this.stage <- srs$stage[i]
  area <- srs$island[i]
  
  # print start 
  print(paste0("Starting ", this.species, " ", this.site, " ", this.stage))
  
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
    filter((latitude_se < 0.05 & longitude_se < 0.125 | #greater allowance for longitude as this can be compressed at poles
             (lon_se_km < 5000 & lat_se_km < 5000)))
  
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
  
  # if no SIC in background, allow
  
  # create binary presence/absence cols
  ars$pa <- 1
  back$pa <- 0
  
  # select key variables
  columns <- c("individual_id", "ed2", "depth", "curr", "pa", "x", "y")
  ars <- ars %>% 
    select(all_of(columns))
  back <- back %>%
    select(all_of(columns))
  
  # reproject ars to epsg:4326
  ars <- ars %>%
    vect(geom = c("x", "y"), crs = "epsg:6932") %>%
    project("epsg:4326") %>%
    as.data.frame(geom = "XY")
  
  # join datasets together
  data <- bind_rows(ars, back)
  
  
  # 3. GAMM
  
  # print gamm initiation
  print("Fitting GAMM")
  
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
  
  # scale covariates 
  data <- data %>%
    mutate(depth = scale(depth)[,1],
           curr = scale(curr)[,1])
  
  # randomise order
  data <- data %>% slice_sample(n = nrow(data))
  
  # run GAMM 
  m1 <- gamm4(pa ~ s(ed2, bs = "ts") + s(depth, bs = "ts") + s(curr, bs = "ts"),
              random = ~(1|individual_id), family = binomial, data = data,
              REML = F)
  
  # get AIC and BIC
  aic_m1 <- AIC(m1$mer)
  bic_m1 <- BIC(m1$mer)
  
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
  
  # get autocorrelation plot
  acf1 <- acf(residuals(m1$gam))
  
  
  # 4. Export
  
  # export model
  saveRDS(m1, paste0("output/gamms/models/", this.species, "/", this.site, " ", this.stage, " gamm.rds"))
  
  # export smooths
  saveRDS(sm, paste0("output/gamms/smooths/", this.species, "/", this.site, " ", this.stage, " smooths.rds"))
  
  # export gamm plot
  ggsave(paste0("output/gamms/plots/", this.species, "/", this.site, " ", this.stage, " gamm.png"), 
         p1, width = 10, height = 6, create.dir = T)
  
  # export autocorrelation plot
  saveRDS(acf1, paste0("output/gamms/acf/", this.species, " ", this.site, " ", this.stage, " acf.rds"))
  
  
  # 5. Model Diagnostics
  
  # print model diagnostics initiation
  print("Calculating delta AIC")
  
  # fit full model without eddies
  m2 <- gamm4(pa ~ s(depth, bs = "ts") + s(curr, bs = "ts"),
              random = ~(1|individual_id), family = binomial, data = data,
              REML = F)
  
  
  # get AIC and BIC
  aic_m2 <- AIC(m2$mer)
  bic_m2 <- BIC(m2$mer)
  
  # make table comparing the two models
  model_comp <- data.frame(model = c("main", "null"),
                           AIC = c(aic_m1, aic_m2),
                           BIC = c(bic_m1, bic_m2))
  
  # export model comparison
  saveRDS(model_comp, paste0("output/gamms/comparison/", this.species, " ", this.site, " ", this.stage, " model comparison.rds"))
  
  # build a separate GAMM excluding 20% of individuals to get CBI metric
  print("Calculating CBI")
  
  # first, randomly select 20% of individuals
  n_inds <- n_distinct(data$individual_id)
  n_exclude <- round(n_inds * 0.2)
  
  # repeat 5 times
  for(i in 1:5){
    print(i)
    test_inds <- sample(unique(data$individual_id), n_exclude)
    
    # remove inds from data
    sub_data <- data %>% filter(!individual_id %in% test_inds)
    
    # get testing data 
    test_data <- data %>% filter(individual_id %in% test_inds)
    
    # fit model to subset of data
    m_sub <- gamm4(pa ~ s(ed2, bs = "ts") + s(depth, bs = "ts") + s(curr, bs = "ts"),
                   random = ~(1|individual_id), family = binomial, data = sub_data,
                   REML = F)
    
    # predict to testing data
    test_data$pred <- predict(m_sub$gam, newdata = test_data, type = "response", exclude = "individual_id")
    
    # calculate CBI
    test_data$pa <- as.factor(test_data$pa)
    test_data$pa <- ordered(test_data$pa, levels = c(1, 0))
    cbi <- boyce_cont(test_data, pa, pred) %>% pull(.estimate)
    
    # join to other cbi
    if(i == 1){
      cbi_df <- data.frame(iteration = i, cbi = cbi)
    } else {
      cbi_df <- bind_rows(cbi_df, data.frame(iteration = i, cbi = cbi))
    }
  }
  
  # get mean CBI across iterations
  cbi <- mean(cbi_df$cbi, na.rm = T)
  
  # export CBI
  saveRDS(cbi_df, paste0("output/gamms/cbi/", this.species, " ", this.site, " ", this.stage, " CBI iterations.rds"))
  saveRDS(cbi, paste0("output/gamms/cbi/", this.species, " ", this.site, " ", this.stage, " CBI.rds"))
  
  # print completion
  print(paste0(this.species, " ", this.site, " ", this.stage, " complete"))
  
}

