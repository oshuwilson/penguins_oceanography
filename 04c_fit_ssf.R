#run step-selection functions

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(terra)
  library(amt)
  library(tidyverse)
  library(mgcv)
  library(gratia)
}


# 1. Setup

#define species site and stage
this.species <- "KIPE"
this.site <- "Crozet"
this.stage <- "incubation"

#read in steps with extracted covariates
steps <- readRDS(paste0("output/ssf/extracted/", this.species, "/", this.site, "_", this.stage, "_steps_extracted.RDS"))

#filter to one individual
ind <- steps %>% filter(individual_id == "Individual_10")

m0 <- ind |> fit_clogit(case_ ~ ssh + curr + eddies + sl_ + cos(ta_) + strata(step_id_))
summary(m0)$coef

m1 <- ind %>%
  fit_issf(case_ ~ ssh + curr + eddies + sl_ + cos(ta_) + strata(step_id_), model = T)
summary(m1)


# 2. Fit GAM

#dummy time variable
steps$times <- 1
ind$times <- 1

#create obs column
steps$obs <- as.numeric(steps$case_)
ind$obs <- as.numeric(ind$case_)

#smooth covariates
m2 <- gam(cbind(times, step_id_) ~ 
            sl_ + log(sl_) + cos(ta_) + 
            s(eddies, k = 5),
          data = ind, 
          method = "REML",
          family = cox.ph,
          weights = obs)
summary(m2)

#hierarchical smooths
m3 <- gam(cbind(times, step_id_) ~ 
            sl_ + log(sl_) + cos(ta_) + 
            s(front_freq, k = 5) + 
            s(front_freq, individual_id, k = 5, bs = "fs"),
          data = steps, 
          method = "REML",
          family = cox.ph,
          weights = obs)
summary(m3)

m4 <- gam(obs ~ sl_ + log(sl_) + cos(ta_) + s(sst, k = 5) + s(ssh, k = 5) + s(mld, k = 5) + s(sal, k = 5) + s(dist_to_colony, k = 5), data = ind, family = binomial)
summary(m4)
