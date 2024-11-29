#hidden markov models in hmmTMB
#are the step lengths for extracted data an issue?? i.e. where NAs were removed???
#split into two scripts and run checks on state assignments before fitting second model

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
  library(miceRanger)
}

# 1. Data preparation

#define variables
this.species <- "ADPE"
this.site <- "Pointe Geologie"
this.stage <- "incubation"

#read in tracks
tracks <- readRDS(file = paste0("output/extractions/", this.species, "/", this.site, "/", this.stage, "_extracted.rds"))
tracks <- tracks %>% 
  filter(pa == "presence")

#read in relevant metadata
meta <- read.csv("~/OneDrive - University of Southampton/Documents/RAATD/RAATD_metadata.csv")
meta <- meta %>%
  filter(abbreviated_name == this.species &
           individual_id %in% unique(tracks$individual_id))

#project tracks to stereographic projection
tracks <- tracks %>% 
  vect(geom = c("x", "y"), crs = "epsg:4326") %>%
  project("epsg:6932") 

#plot to visualise
plot(tracks, pch = ".")

#split tracks into trips using custom function - check plot to see whether buffer distance is appropriate
source("code/functions/trip_split.R")
trax <- trip_split(tracks, meta, buff.dist = 20000)

#format data for hmmTMB
trax <- trax %>%
  rename(ID = trip) %>%
  mutate(ID = as.factor(ID)) %>% 
  as.data.frame()

#covariate-free dataframe
cov_free <- trax %>%
  dplyr::select(c("ID", "x", "y", "date"))

#create hmmTMB object
data <- prepData(cov_free, 
                 type = "UTM", #easting/northing formatting 
                 #covNames = c("depth"), #covariate column names (add later)
)


# 2. Run a covariate-free model to find Par0 values

#establish hidden process model
hid1 <- MarkovChain$new(data = data, n_states = 2)

#define observation distributions
dists <- list(step = "gamma2", angle = "vm")

#visualise step lengths for initial estimates
hist(data$step)

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
  
  #step length mean
  stepMean0 <- runif(2, 
                     min = c(500, 3000),
                     max = c(1500, 8000))
  
  #step length sd
  stepSd0 <- runif(2, 
                   min = c(500, 3000),
                   max = c(1500, 8000))
  
  #turning angle mean
  angleMean0 <- c(0,0)
  
  #turning angle concentration
  angleConc0 <- runif(2, 
                      min = c(0.5, 5),
                      max = c(2, 15))
  
  #combine together
  step <- list(mean = stepMean0, sd = stepSd0)
  angle <- list(mu = angleMean0, kappa = angleConc0)
  par0 <- list(step = step, angle = angle)
  
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
stepMean0 <- as.numeric(top_par[1,])
stepSD0 <- as.numeric(top_par[2,])
angleMean0 <- as.numeric(top_par[3,])
angleConc0 <- as.numeric(top_par[4,])

#combine together
step <- list(mean = stepMean0, sd = stepSD0)
angle <- list(mu = angleMean0, kappa = angleConc0)
par0 <- list(step = step, angle = angle)

# 3. Incorporate covariates to see covariate effects

#vars of interest
vars <- c("depth", "slope", "dshelf", "sst", "mld", "sal", "ssh", "curr", "front_freq", "eddies",
          "dist2ice", "leads", "sic")

#isolate predictors for imputation
preds <- trax %>%
  dplyr::select(all_of(vars))

#impute missing covariate values
if(sum(is.na(preds)) > 0){
  imp <- miceRanger(preds, m = 1)
  preds <- completeData(imp)[[1]]
  
  #re-add imputed data
  trax <- trax %>% 
    dplyr::select(-all_of(vars)) %>%
    bind_cols(preds)
}


#re-prep data
data <- prepData(trax, 
                 type = "UTM", #easting/northing formatting 
                 coordNames = c("x", "y"), #co-ordinate column names
)

#shift covariates up by one index to align with steps and angles
shift_up <- function(col){
  c(col[-1], NA)
}
trax <- trax %>%
  mutate(across(all_of(vars), shift_up))


# 5. Fit HMM with best formula and par0 values

#formula for oceanographic covariates
fixed_formula <- ~ front_freq + s(eddies, k = 5, bs = "ts") + curr +
  s(ID, bs = "re") + s(ID, by = front_freq, bs = "re") + 
  s(ID, by = eddies, bs = "re") + s(ID, by = curr, bs = "re")

#formula for sea ice covariates
ice_formula <- ~ leads + s(dist2ice, k = 5, bs = "ts") + s(sic, k = 5, bs = "ts") + 
  s(ID, bs = "re") + s(ID, by = leads, bs = "re") + 
  s(ID, by = dist2ice, bs = "re") + s(ID, by = sic, bs = "re")

#new hidden process model
hid3 <- MarkovChain$new(data = data, n_states = 2, formula = ice_formula)

#new observation model
obs3 <- Observation$new(data = data, dists = dists, par = par0, n_states = 2)

#create and fit HMM
hmm3 <- HMM$new(obs = obs3, hid = hid3, init = top_hmm)
hmm3$fit(silent = TRUE)

#plot step length and angle distributions
step_lengths <- hmm3$plot_dist("step")
angles <- hmm3$plot_dist("angle")

#plot tracks
hmm_tracks <- hmm3$plot_ts("x", "y") + coord_equal()

#plot covariate effects
curr <- hmm3$plot("tpm", var = "curr")
frfr <- hmm3$plot("tpm", var = "front_freq")
edd <- hmm3$plot("tpm", var = "eddies")
leads <- hmm3$plot("tpm", var = "leads")
sic <- hmm3$plot("tpm", var = "sic")
dist2ice <- hmm3$plot("tpm", var = "dist2ice")


#fit viterbi states to tracks
viterbi <- hmm3$viterbi()
trax <- trax %>%
  mutate(vit = viterbi) %>%
  mutate(state = ifelse(vit == 1, "ARS", "Transit")) %>%
  dplyr::select(-vit)


# 6. Export

longname <- "Adelie Penguins"

#graph formatting
step_lengths <- step_lengths + 
  scale_color_manual("", values = c("darkred", "steelblue4", "black"), labels = c("ARS", "Transit", "Total")) +
  scale_linetype_manual("", values = c(1, 1, 3), labels = c("ARS", "Transit", "Total")) +
  ggtitle(paste0("Step Lengths for ", longname, " at ", this.site, " (", this.stage, ")")) + 
  xlab("Step Length (m)") + ylab("Density")

angles <- angles + 
  scale_color_manual("", values = c("darkred", "steelblue4", "black"), labels = c("ARS", "Transit", "Total")) +
  scale_linetype_manual("", values = c(1, 1, 3), labels = c("ARS", "Transit", "Total")) +
  ggtitle(paste0("Turning Angles for ", longname, " at ", this.site, " (", this.stage, ")")) + 
  xlab("Turning Angle (radians)") + ylab("Density")

hmm_tracks <- hmm_tracks + 
  scale_color_manual("", values = c("red3", "steelblue4"), labels = c("ARS", "Transit")) +
  xlab("Easting (m)") + ylab("Northing (m)") +
  ggtitle(paste0("Tracks for ", longname, " at ", this.site, " (", this.stage, ")"))

curr <- curr + 
  ggtitle(paste0("Relationship between currents and transition probablities for ", longname, " at ", this.site, " (", this.stage, ")")) + 
  xlab("Current (m/s)") + ylab("Transition Probability")

frfr <- frfr + 
  ggtitle(paste0("Relationship between front frequency and transition probablities for ", longname, " at ", this.site, " (", this.stage, ")")) + 
  xlab("Front Frequency") + ylab("Transition Probability")

edd <- edd +
  ggtitle(paste0("Relationship between eddies and transition probablities for ", longname, " at ", this.site, " (", this.stage, ")")) + 
  xlab("Eddies") + ylab("Transition Probability")

leads <- leads + 
  ggtitle(paste0("Relationship between leads and transition probablities for ", longname, " at ", this.site, " (", this.stage, ")")) + 
  xlab("Leads") + ylab("Transition Probability")

sic <- sic +
  ggtitle(paste0("Relationship between sea ice concentration and transition probablities for ", longname, " at ", this.site, " (", this.stage, ")")) + 
  xlab("Sea Ice Concentration") + ylab("Transition Probability")

dist2ice <- dist2ice +
  ggtitle(paste0("Relationship between distance to the ice edge and transition probablities for ", longname, " at ", this.site, " (", this.stage, ")")) + 
  xlab("Distance to Ice (km)") + ylab("Transition Probability")

#list all plots together
plots <- list(step_lengths, angles, hmm_tracks, curr, frfr, edd)
plots <- list(step_lengths, angles, hmm_tracks, leads, sic, dist2ice)

#export plots, tracks with viterbi states, and HMM model
saveRDS(plots, file = paste0("exploration/hmm_outputs/rds/", this.species, "_", this.site, "_", this.stage, "_plots.rds"))
saveRDS(trax, file = paste0("exploration/hmm_outputs/rds/", this.species, "_", this.site, "_", this.stage, "_tracks_with_states.rds"))
saveRDS(hmm3, file = paste0("exploration/hmm_outputs/rds/", this.species, "_", this.site, "_", this.stage, "_hmm.rds"))
