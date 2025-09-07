#---------------------------------------------
# Summarise by species/region
#---------------------------------------------

rm(list = ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

library(tidyverse)

# read in species, site, stage metadata with specialisation info
srs <- read.csv("data/tracks/species_site_stage_v3.csv")

# summarise overall
overall <- srs %>%
  group_by(specialised) %>%
  summarise(n = n())
overall

# summarise by species
by_species <- srs %>% 
  group_by(species, specialised) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = specialised, values_from = n) %>%
  mutate(ratio = yes/(yes + no)) %>%
  ungroup() %>%
  select(species, yes, no, ratio) %>%
  mutate(total = yes + no)
by_species

# summarise by region
by_region <- srs %>% 
  group_by(region, specialised) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = specialised, values_from = n) %>%
  mutate(yes = ifelse(is.na(yes), 0, yes),
         no = ifelse(is.na(no), 0, no)) %>%
  mutate(ratio = yes/(yes + no)) %>%
  ungroup() %>%
  select(region, yes, no, ratio) %>%
  mutate(total = yes + no)
by_region

# summarise by central-place constraints
by_cpf <- srs %>% 
  mutate(cpf = ifelse(cpf %in% c("high", "mid"), "cpf", "not_cpf")) %>%
  group_by(cpf, specialised) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = specialised, values_from = n) %>%
  mutate(yes = ifelse(is.na(yes), 0, yes),
         no = ifelse(is.na(no), 0, no)) %>%
  mutate(ratio = yes/(yes + no)) %>%
  ungroup() %>%
  select(cpf, yes, no, ratio) %>%
  mutate(total = yes + no)
by_cpf


#-------------------------------------------------------------------------
# Chi Squared Tests
#-------------------------------------------------------------------------

# by species
by_species <- by_species %>%
  pivot_longer(cols = c("yes", "no"), names_to = "spec", values_to = "n")

species_cont <- xtabs(n ~ species + spec, data = by_species)
f1 <- fisher.test(species_cont)
c1 <- chisq.test(species_cont)
mosaicplot(species_cont, shade = T)

# by region
by_region <- by_region %>%
  pivot_longer(cols = c("yes", "no"), names_to = "spec", values_to = "n")

region_cont <- xtabs(n ~ region + spec, data = by_region)
c2 <- chisq.test(region_cont)
f2 <- fisher.test(region_cont)
mosaicplot(region_cont, shade = T)


# by cpf
by_cpf <- by_cpf %>%
  pivot_longer(cols = c("yes", "no"), names_to = "spec", values_to = "n")

cpf_cont <- xtabs(n ~ cpf + spec, data = by_cpf)
c3 <- chisq.test(cpf_cont)
f3 <- fisher.test(cpf_cont)
mosaicplot(cpf_cont, shade = T)

# fisher's tests
f1
f2
f3
