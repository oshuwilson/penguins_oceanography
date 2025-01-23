#---------------------------
# ARS Prevalence GAMMs
#---------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(mgcv)
  library(gamm4)
  library(ggside)
  library(patchwork)
  library(gratia)
  library(foreach)
  library(suncalc)
}

#function to fit GAMM using covariates with appropriate unique values
pred_gam <- function(polynyas){
  if(polynyas == FALSE){
    gamm4::gamm4(
      as.formula(
        paste0(
          "bin_state ~ s(", 
          vars %>% paste0(collapse = ", bs = 'ts', k=5) + s("),
          ", bs = 'ts', k=5)")
        ), 
      random = ~(1|individual_id), family=binomial, data=tracks)
  } else {
    gamm4::gamm4(
      as.formula(
        paste0(
          "bin_state ~ s(", 
          vars %>% paste0(collapse = ", bs = 'ts', k=5) + s("),
          ", bs = 'ts', k=5) + polynyas"
        )), 
     random = ~(1|individual_id), family=binomial, data=tracks)
  }
}

# 1. Data Preparation

# define species
this.species <- "MAPE"
longname <- "Macaroni Penguin"

# read in species, region, and stage info
srs <- read.csv("data/tracks/species_site_stage.csv")
srs <- srs %>%
  filter(species == this.species)

# loop over stages
stages <- unique(srs$stage)
for(this.stage in stages){
  
  # loop over regions
  regions <- srs %>%
    filter(stage == this.stage) %>%
    pull(site)
  
  for(this.site in regions){
    try({
    
    # load the tracks with states
    tracks <- readRDS(paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_checked.rds"))
    
    # change state to binomial variable
    tracks <- tracks %>%
      mutate(bin_state = if_else(state == "ARS", 1, 0))
    
    #resample eddy 0s to values between -0.99 and 0.99 - to avoid misleading trends between -1 and 1
    tracks <- tracks %>%
      group_by(row_number()) %>%
      mutate(ed2 = ifelse(eddies == 0, runif(1, -0.99, 0.99), eddies)) %>%
      ungroup()
    
    #read in original tracks to get lat/lons and error info
    original <- readRDS(paste0("output/tracks/", this.species, "/", this.site, " ", this.stage, " tracks.RDS"))
    
    #append latitudes, longitudes, and errors to state tracks
    tracks <- tracks %>% 
      left_join(select(original, individual_id, date, lon, lat, 
                       longitude_se, latitude_se, 
                       lon_se_km, lat_se_km))
    
    #remove tracks with large error
    tracks <- tracks %>%
      filter((latitude_se < 0.133 & longitude_se < 0.133) |
               (lon_se_km < 9000 & lat_se_km < 9000))
    
    #create column in date format for suncalc
    tracks <- tracks %>%
      rename(datetime = date) %>%
      mutate(date = as_date(datetime))
    
    #get sunrise times
    tracks$sunrise <- getSunlightTimes(data = tracks,
                                       keep = c("sunrise"), tz = "UTC") %>%
      pull(sunrise)
    
    #get sunset times
    tracks$sunset <- getSunlightTimes(data = tracks,
                                      keep = c("sunset"), tz = "UTC") %>%
      pull(sunset)
    
    #only keep points between sunrise and sunset
    tracks <- tracks %>%
      filter(datetime >= sunrise & datetime <= sunset |
               is.na(sunrise) & is.na(sunset))
    
    #make polynyas a factor for GAMM
    if("polynyas" %in% names(tracks)){
      tracks <- tracks %>%
        mutate(polynyas = as.factor(polynyas))
    }
    
    #set seed for reproducibility
    set.seed(777)
    
    
    # 2. Fit GAMMs
    
    #list of all possible covariates
    allvars <- c("ed2", "curr", "front_freq", #"depth", "ssh", "mld",
              "dist2ice", "polynyas", "leads", "sic")
    
    #if KIPE or MAPE, remove cryosphere covariates
    if(this.species %in% c("MAPE", "KIPE")){
      allvars <- c("ed2", "curr", "depth", "ssh", "mld", "front_freq")
    }
    
    #find covariates with too little variation for modelling (fewer than 6 unique values)
    lowvar <- tracks %>%
      select(all_of(allvars)) %>% 
      summarise_all(n_distinct) %>%
      pivot_longer(cols = everything(), names_to = "covariate", values_to = "unique") %>%
      filter(unique < 10) %>%
      pull(covariate)
    
    #remove these covariates (excluding polynyas, which is binary by nature)
    if("polynyas" %in% lowvar){
      lowvar <- lowvar[lowvar != "polynyas"]
    }
    vars <- allvars[!allvars %in% lowvar]
    
    #remove covariates with unsuitable variation between individuals (i.e. only one individual shows variation)
    #typically only front_freq or leads
    smallvar <- tracks %>%
      select(front_freq, leads, sic, individual_id) %>%
      group_by(individual_id) %>%
      summarise_all(n_distinct)
    
    ff <- smallvar %>% 
      filter(front_freq == 1) %>% 
      nrow()
    if(ff >= nrow(smallvar) - 1){
      vars <- vars[vars != "front_freq"]
    }
    
    ll <- smallvar %>% 
      filter(leads == 1) %>%
      nrow()
    if(ll >= nrow(smallvar) - 1){
      vars <- vars[vars != "leads"]
    }
    
    ss <- smallvar %>% 
      filter(sic == 1) %>%
      nrow()
    if(ss >= nrow(smallvar) - 1){
      vars <- vars[vars != "sic"]
    }
    
    #if all polynyas are NA or 0, then polynyas = FALSE for pred_gam function
    if(this.species %in% c("ADPE", "CHPE", "EMPE")){
      if(sum(is.na(tracks$polynyas)) == nrow(tracks) |
         nlevels(tracks$polynyas) == 1){
        polynyas = FALSE
      } else{
        polynyas = TRUE
      }
    } else {
      polynyas = FALSE
    }
     
    #now remove polynyas from var list for pred_gam function to work
    vars <- vars[vars != "polynyas"]
    
    #fit GAMM using custom function
    m1 <- pred_gam(polynyas = polynyas)
    
    # get p-values
    m1sum <- as.list(summary(m1$gam))
    pvals <- m1sum$s.pv
    
    # create table of covariate names and p-values
    if(polynyas == TRUE){
      cov_names <- c(vars, "polynyas")
    } else{
      cov_names <- vars
    }
    pscores <- data.frame(cov_names, pvals)
    
    # get smooths
    sm <- smooth_estimates(m1$gam) %>%
      add_confint()
    
    # add constant
    constant <- coef(m1$gam)[1]
    sm <- sm %>%
      mutate(.estimate = .estimate + constant,
             .lower_ci = .lower_ci + constant,
             .upper_ci = .upper_ci + constant)
    
    
    # transform the smooths
    backtrans <- inv_link(m1$gam)
    sm <- sm %>%
      mutate(.estimate = backtrans(.estimate),
             .lower_ci = backtrans(.lower_ci),
             .upper_ci = backtrans(.upper_ci))
    
    # add site and stage to dataframe
    sm <- sm %>%
      mutate(site = this.site,
             stage = this.stage)
    
    # add NAs for any covariates not used in the GAMM
    missing <- allvars[!allvars %in% cov_names]
    for(missingvar in missing){
      sm <- sm %>%
        mutate(!!missingvar := NA)
    }
    
    #return the sm dataframe
    if(this.stage == stages[1] & this.site == regions[1]){
      smooths <- sm
    } else {
      smooths <- rbind(smooths, sm)
    }
    
    #return the p-values as a dataframe
    pscores <- pscores %>%
      mutate(site = this.site,
             stage = this.stage)
    if(this.stage == stages[1] & this.site == regions[1]){
      pvalues <- pscores
    } else {
      pvalues <- rbind(pvalues, pscores)
    }
    })
  }
  
  
  # 3. Plot GAMMs
  
  # front frequency data
  fronts <- smooths %>%
    select(.estimate, .se, .lower_ci, .upper_ci, front_freq, site, stage) %>%
    filter(!is.na(front_freq) &
             stage == this.stage)
  
  # front frequency plot
  frontplot <- ggplot(fronts, aes(x = front_freq, y = .estimate)) +
    geom_line(aes(col = site)) + 
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) + 
    theme_classic() + 
    scale_x_continuous(expand = c(0,0)) +
    scale_color_viridis_d(end = 0.9) +
    scale_fill_viridis_d(end = 0.9) +
    labs(x = "Front Frequency", 
         y = "Probability of ARS",
         fill = "", col = "") +
    theme(legend.position = "top") +
    ggtitle(paste0(longname, "s (", this.stage, ")"))
  
  # eddy data
  eddies <- smooths %>%
    select(.estimate, .se, .lower_ci, .upper_ci, ed2, site, stage) %>%
    filter(!is.na(ed2) &
             stage == this.stage)
  
  # eddy plot
  eddyplot <- ggplot(eddies, aes(x = ed2, y = .estimate)) +
    geom_line(aes(col = site)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) +
    theme_classic() +
    scale_x_continuous(expand = c(0,0)) +
    scale_color_viridis_d(end = 0.9) +
    scale_fill_viridis_d(end = 0.9) +
    labs(x = "Relative Distance to Eddy", 
         y = "Probability of ARS",
         fill = "", col = "") +
    theme(legend.position = "top") +
    ggtitle(paste0(longname, "s (", this.stage, ")"))
  
  
  # current data
  currents <- smooths %>%
    select(.estimate, .se, .lower_ci, .upper_ci, curr, site, stage) %>%
    filter(!is.na(curr) & 
             stage == this.stage)
  
  # current plot
  currplot <- ggplot(currents, aes(x = curr, y = .estimate)) +
    geom_line(aes(col = site)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) +
    theme_classic() +
    scale_x_continuous(expand = c(0,0)) +
    scale_color_viridis_d(end = 0.9) +
    scale_fill_viridis_d(end = 0.9) +
    labs(x = "Current Speed (m/s)", 
         y = "Probability of ARS",
         fill = "", col = "") +
    theme(legend.position = "top") +
    ggtitle(paste0(longname, "s (", this.stage, ")"))
  
  # # depth data
  # depths <- smooths %>%
  #   select(.estimate, .se, .lower_ci, .upper_ci, depth, site, stage) %>%
  #   filter(!is.na(depth) & 
  #            stage == this.stage)
  # 
  # # depth plot
  # depthplot <- ggplot(depths, aes(x = depth, y = .estimate)) +
  #   geom_line(aes(col = site)) +
  #   geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) +
  #   theme_classic() +
  #   scale_x_continuous(expand = c(0,0)) +
  #   scale_color_viridis_d(end = 0.9) +
  #   scale_fill_viridis_d(end = 0.9) +
  #   labs(x = "Depth (m)", 
  #        y = "Probability of ARS",
  #        fill = "", col = "") +
  #   theme(legend.position = "top") +
  #   ggtitle(paste0(longname, "s (", this.stage, ")"))
  # 
  # # ssh data
  # sshs <- smooths %>%
  #   select(.estimate, .se, .lower_ci, .upper_ci, ssh, site, stage) %>%
  #   filter(!is.na(ssh) &
  #            stage == this.stage)
  # 
  # # ssh plot
  # sshplot <- ggplot(sshs, aes(x = ssh, y = .estimate)) +
  #   geom_line(aes(col = site)) +
  #   geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) +
  #   theme_classic() +
  #   scale_x_continuous(expand = c(0,0)) +
  #   scale_color_viridis_d(end = 0.9) +
  #   scale_fill_viridis_d(end = 0.9) +
  #   labs(x = "Sea Surface Height (m)", 
  #        y = "Probability of ARS",
  #        fill = "", col = "") +
  #   theme(legend.position = "top") +
  #   ggtitle(paste0(longname, "s (", this.stage, ")"))
  # 
  # # mld data
  # mlds <- smooths %>%
  #   select(.estimate, .se, .lower_ci, .upper_ci, mld, site, stage) %>%
  #   filter(!is.na(mld) &
  #            stage == this.stage)
  # 
  # # mld plot
  # mldplot <- ggplot(mlds, aes(x = mld, y = .estimate)) +
  #   geom_line(aes(col = site)) +
  #   geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) +
  #   theme_classic() +
  #   scale_x_continuous(expand = c(0,0)) +
  #   scale_color_viridis_d(end = 0.9) +
  #   scale_fill_viridis_d(end = 0.9) +
  #   labs(x = "Mixed Layer Depth (m)", 
  #        y = "Probability of ARS",
  #        fill = "", col = "") +
  #   theme(legend.position = "top") +
  #   ggtitle(paste0(longname, "s (", this.stage, ")"))
  
  # dist2ice data
  if(this.species %in% c("ADPE", "CHPE", "EMPE")){
    if(sum(is.na(smooths$dist2ice)) != nrow(smooths)){
    dist2ices <- smooths %>%
      select(.estimate, .se, .lower_ci, .upper_ci, dist2ice, site, stage) %>%
      filter(!is.na(dist2ice) &
               stage == this.stage)
    
    # dist2ice plot
    dist2iceplot <- ggplot(dist2ices, aes(x = dist2ice, y = .estimate)) +
      geom_line(aes(col = site)) +
      geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) +
      theme_classic() +
      scale_x_continuous(expand = c(0,0)) +
      scale_color_viridis_d(end = 0.9) +
      scale_fill_viridis_d(end = 0.9) +
      labs(x = "Distance to Ice Edge (km)", 
           y = "Probability of ARS",
           fill = "", col = "") +
      theme(legend.position = "top") +
      ggtitle(paste0(longname, "s (", this.stage, ")"))
  }
  
  # lead data
  if(sum(is.na(smooths$leads)) != nrow(smooths)){
    leads <- smooths %>%
      select(.estimate, .se, .lower_ci, .upper_ci, leads, site, stage) %>%
      filter(!is.na(leads) & 
               stage == this.stage)
    
    # lead plot
    leadplot <- ggplot(leads, aes(x = leads, y = .estimate)) +
      geom_line(aes(col = site)) +
      geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) +
      theme_classic() +
      scale_x_continuous(expand = c(0,0)) +
      scale_color_viridis_d(end = 0.9) +
      scale_fill_viridis_d(end = 0.9) +
      labs(x = "Lead Frequency (%)", 
           y = "Probability of ARS",
           fill = "", col = "") +
      theme(legend.position = "top") +
      ggtitle(paste0(longname, "s (", this.stage, ")"))
  }
  
  # polynya plot
  if(sum(is.na(smooths$polynyas)) != nrow(smooths)){
    polynya <- smooths %>%
      select(.estimate, .se, .lower_ci, .upper_ci, polynyas, site, stage) %>%
      filter(!is.na(polynyas) &
               stage == this.stage)
    
    # polynya plot
    polynyaplot <- ggplot(polynya, aes(x = polynyas, y = .estimate)) +
      geom_line(aes(col = site)) +
      geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) +
      theme_classic() +
      scale_x_continuous(expand = c(0,0)) +
      scale_color_viridis_d(end = 0.9) +
      scale_fill_viridis_d(end = 0.9) +
      labs(x = "Polynya Presence", 
           y = "Probability of ARS",
           fill = "", col = "") +
      theme(legend.position = "top") +
      ggtitle(paste0(longname, "s (", this.stage, ")"))
  }
  
  # sic plot
  if(sum(is.na(smooths$sic)) != nrow(smooths)){
    sics <- smooths %>%
      select(.estimate, .se, .lower_ci, .upper_ci, sic, site, stage) %>%
      filter(!is.na(sic) &
               stage == this.stage)
    
    # sic plot
    sicplot <- ggplot(sics, aes(x = sic, y = .estimate)) +
      geom_line(aes(col = site)) +
      geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = site), alpha = 0.2) +
      theme_classic() +
      scale_x_continuous(expand = c(0,0)) +
      scale_color_viridis_d(end = 0.9) +
      scale_fill_viridis_d(end = 0.9) +
      labs(x = "Sea Ice Concentration (%)", 
           y = "Probability of ARS",
           fill = "", col = "") +
      theme(legend.position = "top") +
      ggtitle(paste0(longname, "s (", this.stage, ")"))
  }
  }
  
  #export plots
  ggsave(filename = paste0("output/gamms/plots/", this.species, "/fronts_", this.stage, ".png"),
         plot = frontplot, width = 8, height = 6)
  ggsave(filename = paste0("output/gamms/plots/", this.species, "/eddies_", this.stage, ".png"),
         plot = eddyplot, width = 8, height = 6)
  ggsave(filename = paste0("output/gamms/plots/", this.species, "/currents_", this.stage, ".png"),
         plot = currplot, width = 8, height = 6)
  # ggsave(filename = paste0("output/gamms/plots/", this.species, "/depths_", this.stage, ".png"),
  #        plot = depthplot, width = 8, height = 6)
  # ggsave(filename = paste0("output/gamms/plots/", this.species, "/sshs_", this.stage, ".png"),
  #        plot = sshplot, width = 8, height = 6)
  # ggsave(filename = paste0("output/gamms/plots/", this.species, "/mlds_", this.stage, ".png"),
  #        plot = mldplot, width = 8, height = 6)
  if(exists("dist2iceplot")){
    ggsave(filename = paste0("output/gamms/plots/", this.species, "/dist2ices_", this.stage, ".png"),
           plot = dist2iceplot, width = 8, height = 6)
  }
  if(exists("leadplot")){
    ggsave(filename = paste0("output/gamms/plots/", this.species, "/leads_", this.stage, ".png"),
           plot = leadplot, width = 8, height = 6)
  }
  if(exists("polynyaplot")){
    ggsave(filename = paste0("output/gamms/plots/", this.species, "/polynyas_", this.stage, ".png"),
           plot = polynyaplot, width = 8, height = 6)
  }
  if(exists("sicplot")){
    ggsave(filename = paste0("output/gamms/plots/", this.species, "/sics_", this.stage, ".png"),
           plot = sicplot, width = 8, height = 6)
  }
}

#export smooths
saveRDS(smooths, paste0("output/gamms/smooths/", this.species, "_smooths.rds"))

#export 
saveRDS(pvalues, paste0("output/gamms/pvalues/", this.species, "_pvalues.rds"))
