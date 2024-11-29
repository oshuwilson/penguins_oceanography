#random forest resource selection functions for different suites of covariates
#ADAPT FOR FSLE
#PRODUCE PDPS FOR ICE VARIABLES TOO

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(tidymodels)
  library(tidysdm)
  library(ranger)
  library(miceRanger)
  library(doParallel)
  library(vip)
  library(DALEXtra)
}

#define species, site, and stage
this.species <- "ADPE"
this.site <- "Pointe Geologie"
this.stage <- "incubation"

#read in data with covariates extracted
alldata <- readRDS(paste0("output/extractions/", this.species, "/", this.site, "/", this.stage, "_extracted.rds"))

#define all predictor names
all_predictors <- c("depth", "slope", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr",
                    "eddies", "dist2ice", "front_freq", "leads")


# 1. Format dataframe for model

# base predictor set - Topography and Oceanography
preds <- c("depth", "slope", "dshelf", "sst", "mld", "sal", "ssh", "curr", "eddies", "front_freq")

# if minimum dist2ice is less than 100km, add dist2ice
if(min(alldata$dist2ice, na.rm = TRUE) < 100){
  preds <- c(preds, "dist2ice")
}

# if sea ice concentration ever exceeds 1%, add sea ice concentration
if(max(alldata$sic, na.rm = TRUE) > 0.01){
  preds <- c(preds, "sic")
}

# if sic exceeds 1% and the data are from Austral winter, add leads
if(max(alldata$sic, na.rm = TRUE) > 0.01){
  alldata$month <- month(alldata$date)
  if(any(alldata$month %in% c(5,6,7,8,9,10))){
    preds <- c(preds, "leads")
  }
  alldata <- select(alldata, -month)
}

#select only predictor columns
data <- alldata %>% select(all_of(preds), pa, individual_id)

#remove predictors with more than 10% NA values
data <- data[colSums(is.na(data)) < 0.1*nrow(data)]
predictors <- names(data)
predictors <- subset(predictors, predictors != "pa" & predictors != "individual_id")

#impute missing values
if(sum(is.na(data)) > 0){
  imp <- miceRanger(data, m = 1)
  data <- completeData(imp)[[1]]
}

#set individual ID and presence-absence as factors
data <- data %>% 
  mutate(individual_id = as.factor(individual_id),
         pa = as.factor(pa)) 
  
#relevel presence-absence to ensure presence is the reference level
data$pa <- ordered(data$pa, levels = c("presence", "absence"))


# 2. Hyperparameter Tuning

#enable parallelization
registerDoParallel()
cores <- parallel::detectCores() - 2

#ideal number of folds is 10
v <- 10

#if number of individuals is less than 10, change v to n_ind
n_ind <- length(unique(data$individual_id))
if(n_ind < 10){
  v <- n_ind
}

#create cross-validation folds
folds <- group_vfold_cv(data = data, 
                        group = individual_id, #split training/testing data by individual
                        v = v, #number of folds
                        balance = "observations" #roughly the same number of points in each fold
)

#define formula for modelling
rf_rec <- recipe(pa ~ ., data = data) %>%
  update_role(individual_id, new_role = "ID") #let model know that id is not a predictor

#define random forest settings
rf_mod <- rand_forest() %>%
  set_mode("classification") %>%
  set_engine("ranger", #use ranger package
             importance = "impurity", #gini index for importance
  ) %>%
  set_args(trees = 1000, #1000 trees
           mtry = tune(), #tune mtry
           min_n = 1) #minimum number of samples in a node

#create workflow
rf_wf <- workflow() %>%
  add_model(rf_mod) %>%
  add_recipe(rf_rec)

#define mtry manually
sqrt_pred <- length(predictors) %>% sqrt() %>% floor() #the rounded-down square root of number of predictors
mtry_vals <- c(sqrt_pred - 1, sqrt_pred + 1) #vary mtry by 1 either side of the square root
rf_params <- extract_parameter_set_dials(rf_wf) %>% #update model with these mtry values
  update(mtry = mtry(mtry_vals))

#create tuning grid
rf_grid <- grid_regular(rf_params)

#run models with tuning
rf_tun <- tune_grid(rf_wf,
                    resamples = folds,
                    grid = rf_grid,
                    metrics = sdm_metric_set()) #includes boyce index as a tuning parameter

#get boyce scores for each mtry value
mtry_scores <- collect_metrics(rf_tun, summarize = F) %>%
  filter(.metric == "boyce_cont") 

#extract best model
rf_best <- select_best(rf_tun, metric = "boyce_cont")


# 3. Fit best model to full dataset

#set up model
rf_best_mod <- rand_forest() %>%
  set_engine(engine = "ranger", num.threads = cores, importance = "impurity") %>%
  set_mode("classification") %>%
  set_args(trees = 1000, mtry = rf_best$mtry[1], min_n = 1)

#update workflow
rf_best_wf <- rf_wf %>%
  update_model(rf_best_mod)

#run best model on all data
rf_fit <- rf_best_wf %>%
  fit(data)

#visualise variable importance
rf_fit %>% extract_fit_parsnip() %>% vip(num_features = length(predictors))

#extract variable importance values
var_imp <- rf_fit %>% extract_fit_parsnip() %>% vi()

#get explainer
rf_explainer <- explain_tidymodels(model = rf_fit, 
                                   data = select(data, -pa),
                                   y = as.integer(data$pa),
                                   verbose = T)

#compute partial dependence
pdps <- model_profile(rf_explainer, 
                      variables = c("eddies", "front_freq", "curr", "ssh"),
                      N = 500)

pdps <- model_profile(rf_explainer, 
                      variables = c("eddies", "front_freq", "curr", "ssh", 
                                    "dist2ice", "sic", "leads"),
                      N = 500)

#extract pdp predictive values
pdp_ovr <- as_tibble(pdps$agr_profiles) %>%
  rename(x = `_x_`, yhat = `_yhat_`, var = `_vname_`) %>%
  select(var, x, yhat) %>% 
  mutate(yhat = 1 - yhat) #reverts predictions to presences

#plot pdps
p1 <- ggplot(pdp_ovr, aes(x, yhat)) + 
  geom_line(color = "darkblue", linewidth = 1.2) + 
  facet_wrap(~var, scales = "free_x", nrow = 1) + 
  ylim(0, 1) + 
  theme_bw() +
  ylab("Predicted habitat suitability") + 
  xlab("Predictor values")
p1


# 4. Export key information

# mtry values
mtry_scores <- mtry_scores %>% 
  mutate(best_mtry = ifelse(mtry == rf_best$mtry[1], "Best", "Not")) %>%
  select(mtry, .estimate, best_mtry) %>%
  mutate(species = this.species, site = this.site, stage = this.stage) %>%
  rename(boyce = .estimate)
saveRDS(mtry_scores, file = paste0("output/random_forests/", this.species, "/", this.site, "_", this.stage, "_mtry_scores.rds"))

# partial dependence plots
saveRDS(pdp_ovr, file = paste0("output/random_forests/", this.species, "/", this.site, "_", this.stage, "_pdp.rds"))
ggsave(p1, 
       file = paste0("output/random_forests/", this.species, "/000_", this.site, "_", this.stage, "_pdp.png"),
       width = 10, height = 5, dpi = 300)

# model
saveRDS(rf_fit, file = paste0("output/random_forests/", this.species, "/", this.site, "_", this.stage, "_rf.rds"))

# variable importance scores