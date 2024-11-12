#random forest resource selection functions for different suites of covariates
#ADAPT FOR FSLE
#RUN FOR DIFFERENT COMBINATIONS OF COVARIATES
#remove ice variables for subantarctic species

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(caret)
  library(ranger)
  library(miceRanger)
}

#define species, site, and stage
this.species <- "KIPE"
this.site <- "Macquarie"
this.stage <- "incubation"

#read in data with covariates extracted
alldata <- readRDS(paste0("output/extractions/", this.species, "/", this.site, "/", this.stage, "_extracted.rds"))


# 1. Define model predictor configurations

# Model 1 - Topography and Oceanography
m1 <- c("depth", "slope", "dshelf", "ssh", "curr", "eddies", "front_freq")

#define all predictor names
all_predictors <- c("depth", "slope", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr",
                    "epipelagic_nekton", "upper_meso_nekton", "lower_meso_nekton",
                    "upper_mig_meso_nekton", "lower_mig_meso_nekton", "lower_hmig_meso_nekton",
                    "eddies", "dist2ice", "front_freq", "leads")


#select only predictor columns
data <- alldata %>% select(all_of(all_predictors), pa, individual_id)

#remove predictors with more than 10% NA values
data <- data[colSums(is.na(data)) < 0.1*nrow(data)]
predictors <- names(data)
predictors <- subset(predictors, predictors != "pa" & predictors != "individual_id")

#impute missing values
if(sum(is.na(data)) > 0){
  imp <- miceRanger(data, m = 1)
  data <- completeData(imp)[[1]]
}

#set individual ID as factor
data <- data %>% 
  mutate(individual_id = as.factor(individual_id))

#set number of folds for cross-validation (default 10)
k_number <- 10

#if fewer than 10 individuals in data, set k_number to number of individuals - 1
if(nlevels(data$individual_id) < 10){
  k_number <- nlevels(data$individual_id) - 1
}

#establish folds by ID
folds <- groupKFold(group = data$individual_id, k = k_number)

#find rounded-down square root of number of predictors (for mtry values)
sqrtn <- floor(sqrt(length(predictors)))

#setup parameter grid
param_grid <- expand.grid(mtry = (sqrtn - 1):(sqrtn + 1), 
                          splitrule = "gini", 
                          min.node.size = 1)

#setup cross-validation scheme
cv_scheme <- trainControl(method = "cv",
                          number = length(folds),
                          search = "grid",
                          classProbs = TRUE,
                          sampling = "down",
                          summaryFunction = twoClassSummary,
                          index = folds)

#isolate covariates and pa
X <- data %>% select(all_of(predictors))
Y <- as.factor(data$pa) 

#fit the model
rf <- train(x = X,
            y = Y,
            method = "ranger",
            metric = "ROC",
            trControl = cv_scheme,
            tuneGrid = param_grid,
            importance = "impurity")

#inspect
rf

#variable importance
varImp(rf)

#create partial dependence plots - pmethod can be partdep, apartdep, or plotmo
library(plotmo)
plotmo(rf, pmethod = "plotmo", type = "prob", nresponse = "presence")
?plotmo
