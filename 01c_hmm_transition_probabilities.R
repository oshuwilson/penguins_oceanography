#-------------------------------------------------------------
# Refit HMMs with environmental State Transition Probabilities
#-------------------------------------------------------------

#add random effects???

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(moveHMM)
  library(hmmTMB)
  library(terra)
  library(tidyterra)
  library(tidyverse)
  library(miceRanger)
  library(foreach)
  library(doParallel)
}

#prioritise dplyr select
select <- dplyr::select


# 1. Data preparation

#set species
this.species <- "KIPE"
longname <- "King Penguins"

#read in species, region, and stage info 
srs <- read.csv("data/tracks/species_site_stage.csv")
srs <- srs %>%
  filter(species == this.species)

#set region
regions <- unique(srs$site)
for(this.site in regions){
  
  #set stage
  stages <- srs %>% 
    filter(site == this.site) %>% 
    pull(stage)
  for(this.stage in stages){
    
    #read in HMM quality controlled tracks
    tracks <- readRDS(paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_checked.rds"))
    
    #read in par0 values
    par0 <- readRDS(paste0("output/hmm/hmm_pars/", this.species, "/", this.site, "_", this.stage, "_par0.rds"))
    
    #convert ID to factor
    tracks$ID <- as.factor(tracks$ID)
    
    
    # 2. Incorporate covariates into HMM
    
    #vars of interest - oceanographic only for subantarctic, complete for antarctic
    if(this.species %in% c("KIPE", "MAPE")){ 
      vars <- c("depth", "mld", "ssh", "curr", "front_freq", "eddies", "fsle")
    } else {
      vars <- c("depth", "mld", "ssh", "curr", "front_freq", "eddies", "fsle", "leads", "sic", "dist2ice", "polynyas")
    }
    
    #isolate predictors for imputation
    preds <- tracks %>%
      select(all_of(vars))
    
    #impute missing covariate values
    if(sum(is.na(preds)) > 0){
      imp <- miceRanger(preds, m = 1)
      preds <- completeData(imp)[[1]]
      
      #re-add imputed data
      tracks <- tracks %>% 
        select(-all_of(vars)) %>%
        bind_cols(preds)
    }
    
    #only keep variables of interest and other key columns
    trax <- tracks %>%
      select(ID, individual_id, x, y, date, swim_speed, all_of(vars))
    
    #re-prep data
    data <- prepData(trax, 
                     type = "UTM", #easting/northing formatting 
                     coordNames = c("x", "y"), #co-ordinate column names
    )
    
    #shift covariates up by one index to align with turning angles
    shift_up <- function(col){
      c(col[-1], NA)
    }
    data <- data %>%
      mutate(across(all_of(vars), shift_up))
    
    
    # 3. Fit HMM with best formula and par0 values
    
    #set formula for subantarctic and antarctic species
    if(this.species %in% c("KIPE", "MAPE")){
      formula <- ~ front_freq + s(eddies, k = 5, bs = "ts") + curr + fsle + 
        s(depth, k = 5, bs = "ts") + s(mld, k = 5, bs = "ts") + s(ssh, k = 5, bs = "ts") +
        s(individual_id, bs = "re") + s(individual_id, by = front_freq, bs = "re") + 
        s(individual_id, by = eddies, bs = "re") + s(individual_id, by = curr, bs = "re") +
        s(individual_id, by = fsle, bs = "re") + s(individual_id, by = depth, bs = "re") +
        s(individual_id, by = mld, bs = "re") + s(individual_id, by = ssh, bs = "re")
    } else {
      formula <- ~ front_freq + s(eddies, k = 5, bs = "ts") + curr + fsle + 
        s(depth, k = 5, bs = "ts") + s(mld, k = 5, bs = "ts") + s(ssh, k = 5, bs = "ts") + 
        leads + s(sic, k = 5, bs = "ts") + s(dist2ice, k = 5, bs = "ts") + polynyas +
        s(individual_id, bs = "re")
    }
    
    #new hidden process model
    hid3 <- MarkovChain$new(data = data, n_states = 2, formula = formula, 
                            initial_state = 2)
    
    #set distributions of observed variables
    dists <- list(swim_speed = "gamma2", angle = "vm")
    
    #new observation model
    obs3 <- Observation$new(data = data, dists = dists, par = par0, n_states = 2)
    
    #create and fit HMM
    set.seed(777)
    hmm3 <- HMM$new(obs = obs3, hid = hid3)
    hmm3$fit(silent = T)
    
    
    # 4. Plotting transition probability matrices 
    
    #plot tracks
    hmm_tracks <- hmm3$plot_ts("x", "y") + coord_equal()
    
    #function to plot transition probability matrices
    source("code/functions/tpm.R")
    
    #list individuals
    inds <- unique(data$individual_id)
    
    #plot covariate effects
    curr <- tpm(model = hmm3, var = "curr", inds = inds)
    front <- tpm(model = hmm3, var = "front_freq", inds = inds)
    edd <- tpm(model = hmm3, var = "eddies", inds = inds)
    fsle <- tpm(model = hmm3, var = "fsle", inds = inds)
    depth <- tpm(model = hmm3, var = "depth", inds = inds)
    mld <- tpm(model = hmm3, var = "mld", inds = inds)
    ssh <- tpm(model = hmm3, var = "ssh", inds = inds)
    
    #and for cryosphere variables
    if(this.species %in% c("ADPE", "CHPE", "EMPE")){
      leads <- tpm(model = hmm3, var = "leads", inds = inds)
      sic <- tpm(model = hmm3, var = "sic", inds = inds)
      dist2ice <- tpm(model = hmm3, var = "dist2ice", inds = inds)
      polynya <- tpm(model = hmm3, var = "polynyas", inds = inds)
    }
    
    #polish up plots
    hmm_tracks <- hmm_tracks + 
      scale_color_manual("", values = c("red3", "steelblue4"), labels = c("ARS", "Transit")) +
      xlab("Easting (m)") + ylab("Northing (m)") +
      ggtitle(paste0("Tracks for ", longname, " at ", this.site, " (", this.stage, ")"))
    
    curr <- curr + 
      ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
      xlab("Current (m/s)")
    
    front <- front +
      ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
      xlab("Front Frequency")
    
    edd <- edd +
      ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
      xlab("Eddies")
    
    fsle <- fsle +
      ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
      xlab("FSLE")
    
    depth <- depth +
      ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
      xlab("Depth (m)")
    
    mld <- mld +
      ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
      xlab("Mixed Layer Depth (m)")
    
    ssh <- ssh +
      ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
      xlab("Sea Surface Height (m)")
    
    if(this.species %in% c("ADPE", "CHPE", "EMPE")){
      leads <- leads +
        ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
        xlab("Lead Frequency (%)")
      
      sic <- sic +
        ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
        xlab("Sea Ice Concentration")
      
      dist2ice <- dist2ice +
        ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
        xlab("Distance to Ice Edge (m)")
      
      polynya <- polynya +
        ggtitle(paste0(longname, " at ", this.site, " (", this.stage, ")")) + 
        xlab("Polynyas")
    }
    
    
    # 5. Exporting
    
    #fit viterbi states to tracks
    viterbi <- hmm3$viterbi()
    trax <- trax %>%
      mutate(vit = viterbi) %>%
      mutate(state = ifelse(vit == 1, "ARS", "Transit")) %>%
      dplyr::select(-vit)
    
    #list all plots together
    if(this.species %in% c("KIPE", "MAPE")){
      plots <- list(hmm_tracks, curr, front, edd, fsle, depth, mld, ssh)
    } else {
      plots <- list(hmm_tracks, curr, front, edd, fsle, depth, mld, ssh, leads, sic, dist2ice, polynya)
    }
    
    #export plots, tracks with viterbi states, and HMM model
    saveRDS(plots, file = paste0("output/hmm/tpms/", this.species, "/", this.site, "_", this.stage, "_plots.rds"))
    saveRDS(trax, file = paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_with_states.rds"))
    saveRDS(hmm3, file = paste0("output/hmm/hmm_models/", this.species, "/", this.site, "_", this.stage, "_hmm.rds"))
    
    #print completion
    print(paste0("Completed ", this.species, " at ", this.site, " (", this.stage, ")"))
    
  }
}
