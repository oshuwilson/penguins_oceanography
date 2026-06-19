#----------------------------------
# Assess eddy usage by feeding time 
#----------------------------------

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
srs <- read.csv("data/tracks/species_site_stage_auger.csv")

# keep stages that aren't central-place-foraging
# srs <- srs %>% 
#   filter(stage %in% c("post-breeding", "pre-moult", "post-moult", "fledglings") |
#            species == "KIPE" & stage == "late chick-rearing") # KIPE late chick-rearing is free-roaming

i <- 1

# isolate colony and breeding stage
for(i in 1:nrow(srs)){
  
  this.species <- srs$species[i]
  this.site <- srs$site[i]
  this.stage <- srs$stage[i]
  print(paste0("Starting ", this.species, " ", this.site, " ", this.stage))
  
  # load in hmm checked tracks
  tracks <- readRDS(paste0("output/auger_extractions/", this.species, "/", this.site, "_", this.stage, "_extracted.RDS"))
  
  # if number of distinct individuals is less than 3, skip
  if(n_distinct(tracks$individual_id) < 3){
    next
  }
  
  # 1. Process Data
  
  #remove tracks with large error
  tracks <- tracks %>%
    filter((latitude_se < 0.05 & longitude_se < 0.125 | #greater allowance for longitude as this can be compressed at poles
             (lon_se_km < 5 & lat_se_km < 5)))
  
  #create column in date format for suncalc
  tracks <- tracks %>%
    rename(datetime = date,
           lat = y,
           lon = x) %>%
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
    mutate(ed2 = ifelse(eddies_auger > -1 & eddies_auger < 1, 0, eddies_auger))
  
  # read in background samples
  back <- readRDS(paste0("output/auger_extractions/", this.species, "/", this.site, "_", this.stage, "_background_extracted.RDS"))
  
  # resample non-eddies to 0
  back <- back %>%
    mutate(ed2 = ifelse(eddies_auger > -1 & eddies_auger < 1, 0, eddies_auger))
  
  # create binary presence/absence cols
  ars$pa <- 1
  back$pa <- 0
  
  # select key variables
  ars <- ars %>% 
    select(individual_id, ed2, depth, curr, sic, pa)
  back <- back %>%
    select(individual_id, ed2, depth, curr, sic, pa)
  
  
  # 2. GAMM
  
  # print
  print("Fitting GAMM")
  
  # join datasets together
  data <- bind_rows(ars, back)
  
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
    saveRDS(sm, paste0("output/gamms/smooths/", this.species, "/", this.site, " ", this.stage, " auger smooths.RDS"))
    next
  }
  
  # if sea ice concentration varies in only one or no individuals, don't use as a variable
  smallvar <- ars %>%
    select(sic, individual_id) %>%
    group_by(individual_id) %>%
    summarise_all(n_distinct)
  ss <- smallvar %>% 
    filter(sic == 1) %>%
    nrow()
  if(ss >= nrow(smallvar) - 1){
    use_sic <- FALSE
  } else {
    use_sic <- TRUE
  }
  
  # run GAMM - accounting for variable sea ice use
  if(use_sic == TRUE){
    m1 <- gamm4(pa ~ s(ed2, bs = "ts") + s(depth, bs = "ts") + 
                s(curr, bs = "ts") + s(sic, bs = "ts"),
              random = ~(1|individual_id), family = binomial, data = data)
  } else {
    m1 <- gamm4(pa ~ s(ed2, bs = "ts") + s(depth, bs = "ts") + 
                s(curr, bs = "ts"),
              random = ~(1|individual_id), family = binomial, data = data)
  }
  
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
    labs(x = "Relative Eddy Distance", y = "Odds Ratio")
  
  # get autocorrelation plot
  acf1 <- acf(residuals(m1$gam))
  
  
  # 4. Export
  
  # export model
  saveRDS(m1, paste0("output/auger gamms/models/", this.species, " ", this.site, " ", this.stage, " auger gamm.rds"))
  
  # export smooths
  saveRDS(sm, paste0("output/auger gamms/smooths/", this.species, " ", this.site, " ", this.stage, " auger smooths.rds"))
  
  # export gamm plot
  ggsave(paste0("output/auger gamms/plots/", this.species, " ", this.site, " ", this.stage, " auger gamm.png"), 
         p1, width = 10, height = 6, create.dir = T)
  
  # export autocorrelation plot
  saveRDS(acf1, paste0("output/gamms/acf/", this.species, " ", this.site, " ", this.stage, " acf.rds"))
  
  # 5. Model Diagnostics
  
  # print model diagnostics initiation
  print("Calculating delta AIC")
  
  # fit full model without eddies
  if(use_sic == TRUE){
    m2 <- gamm4(pa ~ s(depth, bs = "ts") + s(curr, bs = "ts") + s(sic, bs = "ts"),
                random = ~(1|individual_id), family = binomial, data = data)
  } else {
    m2 <- gamm4(pa ~ s(depth, bs = "ts") + s(curr, bs = "ts"),
                random = ~(1|individual_id), family = binomial, data = data)
  }
  
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
    if(use_sic == TRUE){
      m_sub <- gamm4(pa ~ s(ed2, bs = "ts") + s(depth, bs = "ts") + s(curr, bs = "ts") + s(sic, bs = "ts"),
                     random = ~(1|individual_id), family = binomial, data = sub_data,
                     REML = F)
    } else {
      m_sub <- gamm4(pa ~ s(ed2, bs = "ts") + s(depth, bs = "ts") + s(curr, bs = "ts"),
                     random = ~(1|individual_id), family = binomial, data = sub_data,
                     REML = F)
    }
    
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