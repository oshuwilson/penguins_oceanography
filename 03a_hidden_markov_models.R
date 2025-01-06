#-------------------------------------------------------------------
#Fit hidden markov models for quality control and initial parameters
#-------------------------------------------------------------------

rm(list = ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(moveHMM)
  library(hmmTMB)
  library(terra)
  library(tidyterra)
  library(foreach)
  library(doParallel)
  library(tidyverse)
  library(CCAMLRGIS)
}

#prioritise dplyr select
select <- dplyr::select

# 1. Data preparation

#define variables
this.species <- "KIPE"

#read in species, region, and stage info 
srs <- read.csv("data/tracks/species_site_stage.csv")
srs <- srs %>%
  filter(species == this.species)

#loop over each study site
regions <- unique(srs$site)
for(this.site in regions){

  #identify stages for this region
  stages <- srs %>% 
    filter(site == this.site) %>% 
    pull(stage)
  
  #loop over each stage
  for(this.stage in stages){

    #read in tracks
    tracks <- readRDS(file = paste0("output/extractions/", this.species, "/", this.site, "_", this.stage, "_extracted.rds"))
    tracks <- tracks %>% 
      filter(pa == "presence")
    
    #read in relevant metadata
    meta <- readRDS("data/metadata.RDS")
    meta <- meta %>%
      filter(abbreviated_name == this.species &
               individual_id %in% unique(tracks$individual_id) &
               device_id %in% unique(tracks$device_id))
    
    #project tracks to stereographic projection
    trax <- tracks %>% 
      vect(geom = c("x", "y"), crs = "epsg:4326") %>%
      project("epsg:6932") 
    
    #plot to visualise
    plot(trax, pch = ".")
    
    #convert back to data frame to calculate swim speeds using eastings/northings
    tracks <- as.data.frame(trax, geom = "XY")
    
    
    # 2. Calculate speed relative to current
    
    #calculate horizontal and vertical relocation speed in m/s
    tracks <- tracks %>%
      group_by(individual_id, device_id) %>%
      mutate(lag_x = lag(x),
             lag_y = lag(y),
             lag_date = lag(date)) %>%
      ungroup() %>%
      mutate(x_diff = lag_x - x,
             y_diff = lag_y - y,
             date_diff = as.numeric(difftime(date, lag_date, units = "secs"))) %>%
      mutate(ux = x_diff/date_diff,
             vx = y_diff/date_diff)
    
    #calculate speed relative to current in m/s
    tracks <- tracks %>%
      mutate(u_diff = abs(ux - uo),
             v_diff = abs(vx - vo)) %>%
      mutate(swim_speed = sqrt((u_diff^2) + (v_diff^2))) %>%
      select(-c(lag_x, lag_y, lag_date, x_diff, y_diff, date_diff, ux, vx))
    
    #visualise swim speeds
    ggplot(tracks, aes(x = swim_speed)) + geom_histogram(boundary = 0) 
    trax <- vect(tracks, geom = c("x", "y"), crs = "epsg:6932")
    ggplot() + geom_spatvector(data = trax, aes(col = swim_speed))
    
    
    # 3. Prepare data for hmmTMB
    
    #split tracks into trips using custom function - check plot to see whether buffer distance is appropriate
    source("code/functions/trip_split.R")
    trax <- trip_split(trax, meta, buff.dist = 10000)
    
    #remove points around colony where localised movements are often misidentified as ARS
    trax <- trax %>%
      filter(dist_to_colony > 10000)
    
    #format data for hmmTMB
    trax <- trax %>%
      rename(ID = trip) %>%
      mutate(ID = as.factor(ID)) %>% 
      as.data.frame()
    
    #covariate-free dataframe
    cov_free <- trax %>%
      select(ID, x, y, date, swim_speed)
    
    #create hmmTMB object
    data <- prepData(cov_free, 
                     type = "UTM") #easting/northing formatting
    
    #convert swim_speed to NA when angle is NA - start and end of trips
    data <- data %>%
      select(-step) %>%
      mutate(swim_speed = case_when(
        is.na(angle) ~ NA,
        TRUE ~ swim_speed))
    
    
    # 4. Hidden Markov Model
    
    #establish hidden process model
    hid1 <- MarkovChain$new(data = data, n_states = 2, 
                            initial_state = 2) #assume initial state is transit
    
    #define observation distributions
    dists <- list(swim_speed = "gamma2", angle = "vm")
    
    #visualise swimming speeds for initial estimates
    hist(data$swim_speed)
    
    #for reproducibility 
    set.seed(777)
    
    #list for Par0 values
    all_pars <- list()
    
    #null variable for loglikelihood results
    all_llk <- NULL
    
    #number of tries with different starting values
    niter <- 25
    
    #set up parallelisation
    cl <- makePSOCKcluster(10)
    registerDoParallel(cl)
    
    #loop through iterations
    all_hmm <- foreach(i = 1:niter, .packages = "hmmTMB") %dopar% {
      
      #this iteration
      iter <- i
      
      #swimming speed mean
      spdMean0 <- runif(2, 
                        min = c(0.2, 1),
                        max = c(0.8, 2))
      
      #swimming speed sd
      spdSd0 <- runif(2, 
                      min = c(0.2, 0.5),
                      max = c(0.5, 1))
      
      #turning angle mean
      angleMean0 <- c(0,0)
      
      #turning angle concentration
      angleConc0 <- runif(2, 
                          min = c(0.5, 5),
                          max = c(2, 15))
      
      #combine together
      swim_speed <- list(mean = spdMean0, sd = spdSd0)
      angle <- list(mu = angleMean0, kappa = angleConc0)
      par0 <- list(swim_speed = swim_speed, angle = angle)
      
      #establish observations
      obs1 <- Observation$new(data = data, dists = dists, par = par0, n_states = 2)
      
      #setup HMM
      hmm1 <- HMM$new(obs = obs1, hid = hid1)
      
      #fit HMM
      hmm1$fit(silent = TRUE)
      
      #return HMM
      return(hmm1)
      
    }
    
    #stop parallelisation
    stopCluster(cl)
    
    #identify the best model by maximum log-likelihood (often multiple best due to numerical stability)
    all_llk <- NULL
    for(i in 1:niter){
      this_hmm <- all_hmm[[i]]
      this_llk <- this_hmm$llk()
      loglk <- data.frame(loglike = this_llk, iteration = i)
      all_llk <- bind_rows(all_llk, loglk)
    }
    
    top_model <- all_llk %>% 
      arrange(desc(loglike)) %>%
      slice(1) %>%
      pull(iteration)
    
    #get corresponding best HMM to use par0s later
    top_hmm <- all_hmm[[top_model]]
    
    #remove all_hmm for processing
    rm(all_hmm)
    
    #extract par0 values from top model
    top_obs <- top_hmm$obs()
    top_par <- top_obs$par()[,,1]
    
    #extract individual components
    spdMean0 <- as.numeric(top_par[1,])
    spdSD0 <- as.numeric(top_par[2,])
    angleMean0 <- as.numeric(top_par[3,])
    angleConc0 <- as.numeric(top_par[4,])
    
    #combine together
    swim_speed <- list(mean = spdMean0, sd = spdSD0)
    angle <- list(mu = angleMean0, kappa = angleConc0)
    par0 <- list(swim_speed = swim_speed, angle = angle)
    
    
    # 4. Extract predicted hidden states for quality control
    
    #get predicted hidden states
    trax$state <- top_hmm$viterbi()
    trax <- trax %>%
      mutate(state = case_when(
        .$state == 1 ~ "ARS",
        .$state == 2 ~ "Transit"
      ))
    
    #convert to terra
    terr <- trax %>%
      as.data.frame() %>%
      vect(geom = c("x", "y"), crs = "epsg:6932")
    
    #list of individual trips
    ids <- unique(terr$ID)
    
    #load coastline
    coast <- readRDS("data/coast_vect.RDS")
    
    #plot for each ID and export to PDF
    pdf(paste0("output/hmm/hmm_checks/", this.species, "/", this.site, "_", this.stage, "_hidden_states.pdf"),
        width = 8, height = 10, pointsize = 16)
    
    #for each trip, plot states
    for(i in ids){
      
      #filter tracks to this individual
      ind_trax <- terr %>%
        filter(ID == i)
      
      #crop coast to track extent
      crop_coast <- crop(coast, ext(ind_trax))
      
      #plot with different colors for states
      p1 <- ggplot() + geom_spatvector(data = crop_coast) + 
        geom_spatvector(data = ind_trax, aes(color = state)) +
        scale_color_manual(values = c("darkred", "steelblue4")) + 
        theme_bw() +
        ggtitle(i)
      
      print(p1)  
    }
    
    #finish pdf export
    dev.off()
    
    
    # 5. Export for next steps
    
    #export Par0 for state-transition probability HMM
    saveRDS(par0, paste0("output/hmm/hmm_pars/", this.species, "/", this.site, "_", this.stage, "_par0.rds"))
    
    #export tracks for state assignment quality control
    saveRDS(trax, paste0("output/hmm/hmm_tracks/", this.species, "/", this.site, "_", this.stage, "_tracks_unchecked.rds"))
    
    #export swimming speed plot
    p1 <- top_hmm$plot_dist("swim_speed")
    p1 <- p1 + ggtitle(paste0(this.species, " ", this.site, " ", this.stage)) +
      xlab("Swimming Speed (m/s)") + ylab("Density") +
      scale_color_manual(labels = c("ARS", "Transit", "Total"), values = c("darkred", "steelblue4", "black")) +
      scale_linetype_manual(labels = c("ARS", "Transit", "Total"), values = c("solid", "solid", "dashed")) +
      labs(col = "State", linetype = "State")
    ggsave(paste0("output/hmm/distributions/", this.species, "/", this.site, "_", this.stage, "_swim_speed.png"), 
           p1, create.dir = T,
           width = 10, height = 6)
    
    #export turning angle plot
    p2 <- top_hmm$plot_dist("angle")
    p2 <- p2 + ggtitle(paste0(this.species, " ", this.site, " ", this.stage)) +
      xlab("Turning Angle (radians)") + ylab("Density") +
      scale_color_manual(labels = c("ARS", "Transit", "Total"), values = c("darkred", "steelblue4", "black")) +
      scale_linetype_manual(labels = c("ARS", "Transit", "Total"), values = c("solid", "solid", "dashed")) +
      labs(col = "State", linetype = "State")
    p2
    ggsave(paste0("output/hmm/distributions/", this.species, "/", this.site, "_", this.stage, "_angle.png"), 
           p2, create.dir = T,
           width = 10, height = 6)
    
    #print completion
    print(paste0("Completed ", this.species, " ", this.site, " ", this.stage))
  }
}
