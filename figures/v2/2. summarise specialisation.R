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

by_species <- srs %>%
  mutate(new = ifelse(new %in% c("specialised", "possible"), "specialised", "other")) %>%
  group_by(species, new) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = new, values_from = n) %>%
  replace_na(list(specialised = 0, other = 0)) %>%
  mutate(ratio = specialised/(specialised + other))
by_species

# summarise by region
by_region <- srs %>% 
  group_by(region, new) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = new, values_from = n) %>%
  replace_na(list(specialised = 0, unspecialised = 0, possible = 0, `not encountered` = 0)) %>%
  mutate(ratio = specialised/(specialised + unspecialised + possible + `not encountered`))
by_region

by_region <- srs %>% 
  mutate(new = ifelse(new %in%  c("specialised", "possible"), "specialised", "other")) %>%
  group_by(region, new) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = new, values_from = n) %>%
  replace_na(list(specialised = 0, other = 0)) %>%
  mutate(ratio = specialised/(specialised + other))
by_region

# summarise by central-place constraints
by_cpf <- srs %>% 
  group_by(cpf, new) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = new, values_from = n) %>%
  replace_na(list(specialised = 0, unspecialised = 0, possible = 0, `not encountered` = 0)) %>%
  mutate(ratio = specialised/(specialised + unspecialised + possible + `not encountered`))
by_cpf


by_cpf <- srs %>% 
  mutate(new = ifelse(new %in% c("specialised", "possible"), "yes", "no")) %>%
  group_by(cpf, new) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = new, values_from = n) %>%
  replace_na(list(yes = 0, no = 0)) %>%
  mutate(ratio = yes/(yes + no))
by_cpf  


#-------------------------------------------------------------------------
# Chi Squared Tests
#-------------------------------------------------------------------------

# by species
by_species <- by_species %>%
  pivot_longer(cols = c("specialised", "other"), names_to = "spec", values_to = "n")

species_cont <- xtabs(n ~ species + spec, data = by_species)
f1 <- fisher.test(species_cont)
c1 <- chisq.test(species_cont)
mosaicplot(species_cont, shade = T)

# by region
by_region <- by_region %>%
  pivot_longer(cols = c("specialised", "other"), names_to = "spec", values_to = "n")

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
