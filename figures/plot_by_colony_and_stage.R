#----------------------------------------------------
# Plot colony/breeding stage trends for each variable
#----------------------------------------------------
# add inset map
# add p-values to plots


rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(ggside)
  library(suncalc)
  library(cowplot)
}

# 1. Setup

# define species
this.species <- "KIPE"
longname <- "King Penguin"

# read in species/site/stage info
srs <- read.csv("data/tracks/species_site_stage_v2.csv")

# filter to this species only
srs <- srs %>% 
  filter(species == this.species)

# read in GAMM smooths and p-values for this species
smooths <- readRDS(paste0("output/GAMMs/smooths/", this.species, "_smooths.RDS"))
pvalues <- readRDS(paste0("output/GAMMs/pvalues/", this.species, "_pvalues.RDS"))

# load in coastline
coast <- readRDS("data/coast_vect.RDS")

# list each region
regions <- unique(srs$site)

# loop over each region
for(this.site in regions){

# identify stages for this region
stages <- srs %>% 
  filter(site == this.site) %>% 
  pull(stage)

# loop over each stage
for(this.stage in stages){


# 2. Process tracks for histograms

# load the tracks with states
tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))

# change state to binomial variable
tracks <- tracks %>%
  mutate(bin_state = if_else(state == "ARS", 1, 0))

#resample eddy 0s to values between -0.99 and 0.99 - to avoid misleading trends between -1 and 1
tracks <- tracks %>%
  group_by(row_number()) %>%
  mutate(ed2 = ifelse(eddies == 0, runif(1, -0.99, 0.99), eddies)) %>%
  ungroup()

#read in original tracks to get lat/lons and error info
area <- srs %>% 
  filter(site == this.site & stage == this.stage) %>% 
  pull(island)
original <- readRDS(paste0("output/tracks/", this.species, "/", area, " ", this.stage, " tracks.RDS"))

#append latitudes, longitudes, and errors to state tracks
tracks <- tracks %>% 
  left_join(select(original, individual_id, date, lon, lat, 
                   longitude_se, latitude_se, 
                   lon_se_km, lat_se_km))

#remove tracks with large error
tracks <- tracks %>%
  filter((latitude_se < 0.133 & longitude_se < 0.133) |
           (lon_se_km < 9000 & lat_se_km < 9000))

#create column in date format for suncalc
tracks <- tracks %>%
  rename(datetime = date) %>%
  mutate(date = as_date(datetime))

#get sunrise times
tracks$dawn <- getSunlightTimes(data = tracks,
                                   keep = c("dawn"), tz = "UTC") %>%
  pull(dawn)

#get sunset times
tracks$dusk <- getSunlightTimes(data = tracks,
                                  keep = c("dusk"), tz = "UTC") %>%
  pull(dusk)

#only keep points between sunrise and sunset
tracks <- tracks %>%
  filter(datetime >= dawn & datetime <= dusk |
           is.na(dawn) & is.na(dusk))

# get p-values for this site and stage
pvals <- pvalues %>%
  filter(site == this.site & stage == this.stage)
edp <- pvals %>%
  filter(cov_names == "ed2") %>%
  pull(pvals)


# # 3. Plot tracks
# 
# # convert original tracks to terra
# original <- original %>%
#   vect(geom = c("lon", "lat"), crs = "epsg:4326") %>%
#   project("epsg:6932")
# 
# # crop coast to track extent
# crop_coast <- crop(coast, ext(original))
# 
# # project coast and tracks to lat/lon
# original <- project(original, "epsg:4326")
# crop_coast <- project(crop_coast, "epsg:4326")
# 
# # find the average yday and year of tracks
# dates <- unique(as_date(original$date))
# years <- unique(year(original$date))
# 
# # get current data (vo and uo) for all years
# for(year in years){
#   try(uo <- rast(paste0("E:/Satellite_Data/daily/uo/uo_", year, ".nc")))
#   try(vo <- rast(paste0("E:/Satellite_Data/daily/vo/vo_", year, ".nc")))
#   
#   # if no vo or uo for this year, skip
#   if(!exists("uo") | !exists("vo")){
#     years <- years[years != year]
#     next
#   }
#   
#   # limit to dates in tracks
#   uo <- uo[[time(uo) %in% dates]]
#   vo <- vo[[time(vo) %in% dates]]
#   
#   # crop to extent of tracks 
#   uo <- crop(uo, ext(original) + c(2, 2, 2, 2))
#   vo <- crop(vo, ext(original) + c(2, 2, 2, 2))
#   
#   # calculate current layer
#   curr <- sqrt(uo^2 + vo^2)
# 
#   # join to other years
#   if(year == years[1]){
#     current <- curr
#   } else {
#     current <- c(current, curr)
#   }
#   
#   #reset
#   rm(uo, vo, curr)
# }
# 
# # calculate average current speed
# current <- mean(current)
# 
# # plot coast and tracks
# p1 <- ggplot() +
#   geom_spatraster(data = current) +
#   geom_spatvector(data = crop_coast, fill = "grey") + 
#   geom_spatvector(data = original, size = 0.5, col = "black") +
#   theme_classic() +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   labs(fill = "Mean Current Velocity (m/s)") +
#   theme(legend.position = "top") +
#   scale_fill_continuous(guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
# p1


# 4. Plot Relative Eddy Distance

# get eddy distance
eddies <- smooths %>%
  filter(!is.na(ed2) &
           site == this.site &
           stage == this.stage)

# create eddie background shading values
ed_vals <- data.frame(type = as.factor(c("anticyclonic core", "anticyclonic periphery", "not eddy",
                               "cyclonic periphery", "cyclonic core")),
                      x = c(-2.5, -1.5, 0, 1.5, 2.5),
                      xmin = c(-3, -2, -1, 1, 2),
                      xmax = c(-2, -1, 1, 2, 3), 
                      ymin = 0, ymax = 1)

# plot eddy GAMM output
p3 <- ggplot(eddies) + 
  geom_tile(data = ed_vals, alpha = 0.4, height = Inf, width = c(1, 1, 2, 1, 1),
            aes(x = x, y = 0.5, fill = type)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = ed2),
              alpha = 0.2) +
  geom_line(aes(x = ed2, y = .estimate), lwd = 0.5) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) + 
  ylab("Prevalance of ARS state") +
  xlab("Relative Eddy Distance")

# add histogram
p3 <- p3 + geom_xsidehistogram(data = tracks, aes(x=ed2), 
                               fill = "#230493", alpha = 0.2,
                               boundary = -3, bins = 60) + 
  theme_ggside_void() + 
  theme(ggside.panel.scale.x = .3,
        legend.title = element_blank(),
        legend.position = "bottom")

# add eddy shading
p3 <- p3 + scale_fill_manual(values = c("red4", "red3",
                                        "steelblue4", "steelblue3",
                                        "white"),
                             breaks = c("anticyclonic core", "anticyclonic periphery",
                                        "cyclonic core", "cyclonic periphery", "not eddy")) +
  scale_xsidey_continuous(expand = c(0,0))

# add title
p3 <- p3 + ggtitle(paste0(longname, "s at ", this.site, " (", this.stage, ")"))

# export
ggsave(filename = paste0("output/imagery/colony plots/", this.species, "/", this.site, " ", this.stage, " eddy.png"),
       plot = p3, width = 8, height = 8, units = "in", dpi = 300)


# # remove legend from p3
# p3 <- p3 + guides(fill = "none")
# 
# # get legend
# legend <- readRDS("output/imagery/colony plots/legend.RDS")
# 
# 
# 
# # 6. layout all together
# 
# # oceanography plots
# plots <- plot_grid(p2, p3, legend, nrow = 3, rel_heights = c(4, 4, 1))
# plots
# 
# #all plots
# layout <- plot_grid(p1, plots, ncol = 2)
# layout
# 
# #title
# title <- ggdraw() + 
#   draw_label(
#     paste0(longname, "s at ", this.site, " (", this.stage, ")"),
#     fontface = 'bold',
#     x = 0,
#     hjust = 0
#   ) +
#   theme(
#     # add margin on the left of the drawing canvas,
#     # so title is aligned with left edge of first plot
#     plot.margin = margin(0, 0, 0, 200)
#   )
# 
# #all together
# final <- plot_grid(title, layout, nrow = 2, rel_heights = c(1, 20))
# final
# 
# # export 
# ggsave(paste0("output/imagery/colony plots/", this.species, "/", this.site, "_", this.stage, ".svg"), final, 
#        width = 9, height = 8, units = "in", dpi = 300,
#        create.dir = T)

# print completion
print(paste0("Completed ", this.site, " ", this.stage))

}
}

