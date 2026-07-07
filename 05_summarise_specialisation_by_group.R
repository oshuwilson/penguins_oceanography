#-------------------------------------------------------------------------------
# Summarise results by species, region, and breeding constraints
#-------------------------------------------------------------------------------

rm(list = ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

library(tidyverse)

# read in species, site, stage metadata with specialisation info
srs <- read.csv("data/tracks/species_site_stage_v3.csv")

# summarise overall stats
overall <- srs %>%
  group_by(new) %>%
  summarise(n = n())
overall

# summarise by species
by_species <- srs %>%
  group_by(species, new) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = new, values_from = n) %>%
  replace_na(list(specialised = 0, unspecialised = 0, possible = 0, `not encountered` = 0)) %>%
  mutate(ratio = specialised/(specialised + unspecialised + possible + `not encountered`))
by_species

# summarise by region
by_region <- srs %>% 
  group_by(region, new) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = new, values_from = n) %>%
  replace_na(list(specialised = 0, unspecialised = 0, possible = 0, `not encountered` = 0)) %>%
  mutate(ratio = specialised/(specialised + unspecialised + possible + `not encountered`))
by_region

# summarise by central-place constraints
by_cpf <- srs %>% 
  group_by(cpf, new) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = new, values_from = n) %>%
  replace_na(list(specialised = 0, unspecialised = 0, possible = 0, `not encountered` = 0)) %>%
  mutate(ratio = specialised/(specialised + unspecialised + possible + `not encountered`))
by_cpf


#-------------------------------------------------------------------------
# Fisher's Exact Tests
#-------------------------------------------------------------------------

# by species
by_species <- by_species %>%
  mutate(other = sum(unspecialised, possible, `not encountered`)) %>%
  pivot_longer(cols = c("specialised", "other"), names_to = "spec", values_to = "n")

# contingency table
species_cont <- xtabs(n ~ species + spec, data = by_species)
mosaicplot(species_cont, shade = T)

# Fisher's exact test
f1 <- fisher.test(species_cont)
f1

# by region
by_region <- by_region %>%
  mutate(other = sum(unspecialised, possible, `not encountered`)) %>%
  pivot_longer(cols = c("specialised", "other"), names_to = "spec", values_to = "n")

# contingency table
region_cont <- xtabs(n ~ region + spec, data = by_region)
mosaicplot(region_cont, shade = T)

# Fisher's exact test
f2 <- fisher.test(region_cont)
f2

# by central-place constraints
by_cpf <- by_cpf %>%
  mutate(other = sum(unspecialised, possible, `not encountered`)) %>%
  pivot_longer(cols = c("specialised", "other"), names_to = "specialisation", values_to = "n")

# contingency table
cpf_cont <- xtabs(n ~ cpf + specialisation, data = by_cpf)
mosaicplot(cpf_cont, shade = T)

# Fisher's exact test
f3 <- fisher.test(cpf_cont)
f3
