#------------------------------------------------------
# Plot Background Availability and Eddy Usage by Colony
#------------------------------------------------------

# only for case studies modified following Matthis Auger and JB Sallee's discussion

setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(gamm4)
  library(gratia)
  library(cowplot)
}

rm(list=ls())

# read in species site stage info to loop over
srs <- read.csv("data/tracks/species_site_stage_v2.csv")

# read in updated stage versions
srs2 <- read.csv("data/tracks/species_site_stage_v3.csv")

# read in list of examples using new dataset
auger <- read.csv("data/tracks/species_site_stage_auger.csv") %>%
  mutate(case_study = paste(species, site, stage, sep = " "))

# list new smooth files
new.smooths <- list.files(path = "output/gamms/smooths_pack_ice")

# for each smooth file
this.smooth <- new.smooths[28]

# print smooth
print(this.smooth)

# get species, site, stage, and new.stage 
this.species <- "EMPE"
this.site <- "Taylor Glacier"
this.stage <- "fledglings"
new.stage <- "fledglings"

# does this use the updated eddy product
if(paste(this.species, this.site, this.stage, sep = " ") %in% auger$case_study){
  new_eddies <- TRUE
} else {
  new_eddies <- FALSE
}

#-------------------------------------------------------------------------------
# 1. Background Availability Plot
#-------------------------------------------------------------------------------

# load in hmm checked tracks
if(new_eddies == FALSE){
  tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))
} else {
  tracks <- readRDS(paste0("output/auger_extractions/", this.species, "/", this.site, "_", this.stage, "_extracted.rds"))
}

# if fewer than 3 individuals, skip
if(length(unique(tracks$individual_id)) < 3){
  rm(list=ls())
  print("SKIP")
}


# read in background samples
if(new_eddies == FALSE){
  back <- readRDS(paste0("output/background/", this.species, "/", this.site, " ", this.stage, " background.rds"))
} else {
  back <- readRDS(paste0("output/auger_extractions/", this.species, "/", this.site, "_", this.stage, "_background_extracted.rds"))
}

# load in metadata for colony location
meta <- readRDS("data/metadata.rds")

# filter metadata to this colony
meta <- meta %>% filter(abbreviated_name == this.species & individual_id %in% unique(tracks$individual_id))

# get colony location
colony <- meta %>%
  slice(1) %>%
  select(lon = deployment_decimal_longitude, lat = deployment_decimal_latitude) %>%
  vect(geom = c("lon", "lat"), crs = "epsg:4326") 

# convert background samples and tracks to terra
back <- back %>%
  vect(geom = c("x", "y"), crs = "epsg:4326")
if(new_eddies == F){
tracks <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")
} else {
tracks <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:4326")
}

# calculate maximum distance to colony
back$distances <- distance(back, colony)
tracks$dist_to_colony <- distance(tracks, colony)

# resample eddy values between -1 and 1 to 0
if(new_eddies == F){
  back <- back %>%
    mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))
} else {
  back <- back %>%
    mutate(ed2 = ifelse(eddies_auger > -1 & eddies_auger < 1, 0, eddies))
}

# classify eddies into anticyclonic, cyclonic, or non-eddy
back <- back %>%
  mutate(ed3 = case_when(
    ed2 < 0 ~ "anticyclonic",
    ed2 > 0 ~ "cyclonic",
    ed2 == 0 ~ "non-eddy"
  )) 

# convert distance to km
back <- back %>%
  mutate(distances = distances / 1000)
tracks <- tracks %>%
  mutate(distances = dist_to_colony/ 1000)

# remove NA values for ed2
back <- back %>%
  filter(!is.na(ed2))

# plot background eddy availability
p1 <- ggplot(back, aes(x = distances)) + 
  geom_density(aes(group = ed3, fill = ed3), 
               adjust = 1.5, position = "fill", color = NA) + 
  theme_classic() + 
  labs(x = "Distance to Colony (km)", y = "Relative Area") +
  scale_fill_manual(values = c("darkred", "steelblue4", "lightgrey"),
                    name = "",
                    labels = c("Anticyclonic", "Cyclonic", "Non-Eddy")) +
  stat_density(data = tracks, aes(x = distances, y = after_stat(scaled), 
                                  linetype = "Track Locations"),
               lwd = 1, geom = "line") +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(override.aes = list(linetype = c(0, 0, 0)))) +
  scale_linetype_manual(values = "dashed", name = NULL) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.position = "top")

# if no cyclonic eddies, plot with only two colours
if(max(back$ed2) == 0 & min(back$ed2) < 0){
  p1 <- ggplot(back, aes(x = distances)) + 
    geom_density(aes(group = ed3, fill = ed3), 
                 adjust = 1.5, position = "fill", color = NA) + 
    theme_classic() + 
    labs(x = "Distance to Colony (km)", y = "Relative Area") + 
    scale_fill_manual(values = c("darkred", "lightgrey"),
                      name = "",
                      labels = c("Anticyclonic", "Non-Eddy")) +
    stat_density(data = tracks, aes(x = distances, y = after_stat(scaled), 
                                    linetype = "Track Locations"),
                 lwd = 1, geom = "line") +
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +
    guides(fill = guide_legend(override.aes = list(linetype = c(0, 0)))) +
    scale_linetype_manual(values = "dashed", name = NULL) +
    theme(legend.key.size = unit(1, 'cm'),
          legend.position = "top")
}

# if no anticyclonic eddies, plot with two colours
if(max(back$ed2 > 0) & min(back$ed2) == 0){
  p1 <- ggplot(back, aes(x = distances)) + 
    geom_density(aes(group = ed3, fill = ed3), 
                 adjust = 1.5, position = "fill", color = NA) + 
    theme_classic() + 
    labs(x = "Distance to Colony (km)", y = "Relative Area") + 
    scale_fill_manual(values = c("steelblue4", "lightgrey"),
                      name = "",
                      labels = c("Cyclonic", "Non-Eddy")) +
    stat_density(data = tracks, aes(x = distances, y = after_stat(scaled), 
                                    linetype = "Track Locations"),
                 lwd = 1, geom = "line") +
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +
    guides(fill = guide_legend(override.aes = list(linetype = c(0, 0)))) +
    scale_linetype_manual(values = "dashed", name = NULL) +
    theme(legend.key.size = unit(1, 'cm'),
          legend.position = "top")
}

# if all values are non-eddy print grey plot
if(max(back$ed2) == 0 & min(back$ed2) == 0){
  p1 <- ggplot(back, aes(x = distances)) + 
    geom_density(aes(group = ed3, fill = ed3), 
                 adjust = 1.5, position = "fill", color = NA) + 
    theme_classic() +
    labs(x = "Distance to Colony (km)", y = "Relative Area") + 
    scale_fill_manual(values = c("lightgrey"),
                      name = "",
                      labels = c("Non-Eddy")) +
    stat_density(data = tracks, aes(x = distances, y = after_stat(scaled), 
                                    linetype = "Track Locations"),
                 lwd = 1, geom = "line") +
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +
    guides(fill = guide_legend(override.aes = list(linetype = c(0)))) +
    scale_linetype_manual(values = "dashed", name = NULL) +
    theme(legend.key.size = unit(1, 'cm'),
          legend.position = "top")
}

p1

#-------------------------------------------------------------------------------
# 2. Anticyclone and Cyclone Odds Ratios
#-------------------------------------------------------------------------------

# cleanup
rm(list=setdiff(ls(), c("p1", "i", "srs", "srs2", "this.species", "this.site", "this.stage", "new.stage", "new_eddies", "this.smooth")))

# load in smooths
if(new_eddies == F){
  sm <- readRDS(paste0("output/gamms/smooths_pack_ice/", this.species, " ", this.site, " ", this.stage, " smooths.rds"))
} else {
  sm <- readRDS(paste0("output/gamms/smooths_pack_ice/", this.species, " ", this.site, " ", this.stage, " auger smooths.rds"))
}

# average odds ratios over .5 intervals
sm2 <- sm %>%
  mutate(ed2 = cut(ed2, breaks = c(-3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.125, 
                                   0.125, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3))) %>%
  group_by(ed2) %>%
  summarise(.estimate = mean(.estimate),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci))

# remove two categories with no associated eddy region
sm2 <- sm2 %>% 
  filter(!ed2 %in% c("(-1,-0.125]", "(0.125,1]"))

# get upper and lower limits from cuts
sm2 <- sm2 %>%
  mutate(x_tmp = str_sub(ed2, 2, -2)) %>%
  separate(x_tmp, c("min", "max"), sep = ",") %>%
  mutate_at(c("min", "max"), as.numeric)

# rename odds ratio
sm2 <- sm2 %>%
  rename(OR = .estimate)

# if lower and upper ci envelope 1, set odds ratio to 1
sm2 <- sm2 %>%
  mutate(OR = ifelse(.lower_ci < 1 & .upper_ci > 1, 1, OR))

# read in blank odds ratios if some missing
blanksm <- readRDS("output/imagery/colony odds ratios/sm blank.RDS")
sm2 <- bind_rows(sm2, blanksm) %>%
  distinct(ed2, .keep_all = T)

# create cyclonic dataframe
cyclones <- sm2 %>% 
  filter(min >= 1) %>%
  arrange(desc(min))

# create factor for plotting
cyclones <- cyclones %>%
  mutate(band = factor(1:nrow(cyclones)))

# truncate odds ratios above 3
cyclones <- cyclones %>%
  mutate(OR = ifelse(OR > 3, 3, OR))

# create anticyclonic dataframe
anticyclones <- sm2 %>% 
  filter(max <= -1) %>%
  arrange(min)

# create factor for plotting
anticyclones <- anticyclones %>%
  mutate(band = factor(1:nrow(anticyclones)))

# truncate odds ratios above 3
anticyclones <- anticyclones %>%
  mutate(OR = ifelse(OR > 3, 3, OR))

# plot cyclones
p2 <- ggplot(cyclones, aes(x = band, fill = as.numeric(OR))) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 8.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 1,
                       name = "Odds Ratio", limits = c(0, 3),
                       labels = c(0, 1, 2, "3+"),
                       na.value = "grey") +
  ggtitle("Cyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5)) 
p2

# plot anticyclones
p3 <- ggplot(anticyclones, aes(x = band, fill = as.numeric(OR))) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 8.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 1,
                       name = "Odds Ratio", limits = c(0, 3),
                       labels = c(0, 1, 2, "3+"),
                       na.value = "grey") +
  ggtitle("Anticyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5))
p3

#-------------------------------------------------------------------------------
# 3. Plot all together
#-------------------------------------------------------------------------------

# cleanup
rm(list=setdiff(ls(), c("p1", "p2", "p3", "i", "srs", "srs2", "this.species", "this.site", "this.stage", "new.stage",
                        "new_eddies", "this.smooth")))

# get legend for eddy plots
legend <- get_legend(p3)

# remove legend from eddy plots
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

# plot eddy plots side by side
odds_ratios <- plot_grid(p3, legend, p2, ncol = 3, rel_widths = c(1, 0.05, 1)) 
odds_ratios

# plot background availability above
plots <- plot_grid(p1, odds_ratios, ncol = 1)
plots

# add title
title <- ggdraw() + 
  draw_label(
    paste0(this.species, " ", this.site, " (", new.stage, ")"),
    fontface = 'bold',
    size = 20,
    x = 0,
    hjust = -0.01
  ) 

# all together
final <- plot_grid(title, plots, nrow = 2, rel_heights = c(1, 10), align = "v")
print(final)

# export plots
ggsave(paste0("output/imagery/colony odds ratios pack ice/", this.species, "/", this.site, " ", new.stage, " odds ratios.png"), final, 
       width = 7.8, height = 11, dpi = 300, create.dir = T)

