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
  library(CCAMLRGIS)
}

# 1. Data preparation

#define variables
this.species <- "KIPE"
this.site <- "Marion"
this.stage <- "chick-rearing"

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
trax <- trip_split(tracks, meta, buff.dist = 10000)

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
hid1 <- MarkovChain$new(data = data, n_states = 2, 
                        initial_state = 2) #assume initial state is transit

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
coast <- vect(load_Coastline())

#plot for each ID and export to PDF
pdf(paste0("output/hmm/hmm_checks/", this.species, "_", this.site, "_", this.stage, "_hidden_states.pdf"),
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
saveRDS(par0, paste0("output/hmm/hmm_pars/", this.species, "_", this.site, "_", this.stage, "_par0.rds"))

#export tracks for state assignment quality control
saveRDS(trax, paste0("output/hmm/hmm_tracks/", this.species, "_", this.site, "_", this.stage, "_tracks_unchecked.rds"))

#ADD STEPS AND ANGLES PLOTS
