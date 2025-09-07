#--------------------------------------
# Schematic Plots
#--------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)
library(terra)
library(tidyterra)
library(suncalc)
library(gratia)
library(cowplot)

#-------------------------------------------------------------------------
# 1. Plot of all tracks
#-------------------------------------------------------------------------

# read in tracks for incubating macaroni penguins from Fairy Point
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/MAPE/Fairy Point, Bird Island incubation tracks checked.rds")

# convert to terra
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# create lines for each individual
inds <- unique(trax$individual_id)
for(ind in inds){
  
  ind_trax <- trax %>%
    filter(individual_id == ind)
  start_date <- min(ind_trax$date)
  ind_trax <- as.lines(ind_trax)
  ind_trax$date <- start_date
  ind_trax$individual_id <- ind
  if(ind == inds[1]){
    trax_lines <- ind_trax
  } else {
    trax_lines <- c(trax_lines, ind_trax)
  }
}
trax_lines <- vect(trax_lines)

# get coastline from rnaturalearth
coast <- rnaturalearth::ne_countries(scale = 10, returnclass = "sv")
coast <- crop(coast, ext(trax) + c(0.5, 0.5, 0.5, 0.5))
plot(coast)

# plot
p1 <- ggplot() +
  geom_spatvector(data = coast, col = "white", fill = "white") +
  geom_spatvector(data = trax_lines, aes(col = individual_id)) +
  scale_color_viridis_d() +
  guides(col = "none") +
  theme_void() +
  theme(
    panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
p1
ggsave("output/schematic/tracks.png", p1, bg = "transparent")


#-------------------------------------------------------------------
# 2. Plot hidden states for one track
#-------------------------------------------------------------------

# plot example track
p2 <- ggplot() +
  geom_spatvector(data = coast, col = "grey40", fill = "grey40") +
  geom_spatvector(data = trax %>% 
                    filter(individual_id == "MAPE-dtsetBirdLife751-X14-RAATD_2" &
                             date < as_date("1999-12-30")), 
                  aes(col = state)) +
  scale_color_manual(values = c("red3", "steelblue4"), name = "State") +
  guides(color = "none") +
  theme_void() +
  theme(
    panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent', color = NA), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = NA) #transparent legend panel
  )
p2

# save
ggsave("output/schematic/hmm_states.png", p2, bg = "transparent")


#--------------------------------------------------------------------
# 3. Background Sampling
#--------------------------------------------------------------------

# read in background samples
back <- readRDS("output/background/MAPE/Fairy Point, Bird Island incubation background.rds")

# convert to terra
back <- vect(back, geom = c("x", "y"), crs = "epsg:4326")

# plot 
p3 <- ggplot() +
  geom_spatvector(data = coast, col = "grey40", fill = "grey40") +
  geom_spatvector(data = back, col = "steelblue4", size = 0.5) +
  guides(color = "none") +
  theme_void() +
  theme(
    panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent', color = NA), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = NA) #transparent legend panel
  )
p3

# save
ggsave("output/schematic/bg.png", p3, bg = "transparent")


#------------------------------------------------------------------
# 4. Covariates
#------------------------------------------------------------------

# read in eddies
eddies <- rast("E:/Satellite_Data/daily/eddies/eddies_1999.nc")

# limit to dates in tracking data
dates <- unique(as_date(tracks$date))
eds <- eddies[[time(eddies) %in% dates]]

# crop
eds <- crop(eds, ext(trax) + c(0.5, 0.5, 0.5, 0.5))

# pick example
eds <- eds[[time(eds) == as_date("1999-12-14")]]
time(eds)

# plot
p4 <- ggplot() +
  geom_spatraster(data = eds) +
  geom_spatvector(data = coast, fill = "white") +
  scale_fill_gradient2() +
  guides(fill = "none") +
  theme_void() +
  theme(
    panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent', color = NA), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = NA) #transparent legend panel
  )

# save
ggsave("output/schematic/eddies.png", p4, bg = "transparent")

rm(eddies)

# read in depth
depth <- rast("E:/Satellite_Data/static/depth/depth.nc")

# crop
depth <- crop(depth, ext(trax) + c(0.5, 0.5, 0.5, 0.5))

p5 <- ggplot() +
  geom_spatraster(data = depth) +
  geom_spatvector(data = coast, fill = "white") +
  scale_fill_viridis_c(direction = -1, na.value = "#fde725ff") +
  guides(fill = "none") +
  theme_void() +
  theme(
    panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent', color = NA), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = NA) #transparent legend panel
  )
p5

# save
ggsave("output/schematic/depth.png", p5, bg = "transparent")


# read in current
uo <- rast("E:/Satellite_Data/daily/uo/uo_1999.nc")
vo <- rast("E:/Satellite_Data/daily/vo/vo_1999.nc")

# limit time
uo <- uo[[time(uo) == as_date("1999-12-14")]]
vo <- vo[[time(vo) == as_date("1999-12-14")]]

# crop
e <- ext(trax) + c(0.5, 0.5, 0.5, 0.5)
uo <- crop(uo, e)
vo <- crop(vo, e)

curr <- sqrt(uo^2 + vo^2)
plot(curr)

p6 <- ggplot() +
  geom_spatraster(data = curr) +
  geom_spatvector(data = coast, fill = "white") +
  scale_fill_viridis_c(na.value = "#440154ff") +
  guides(fill = "none") +
  theme_void() +
  theme(
    panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent', color = NA), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = NA) #transparent legend panel
  )
p6
ggsave("output/schematic/curr.png", p6, bg = "transparent")


#-----------------------------------------------------------
# 5. GAMM plots
#-----------------------------------------------------------

# load in GAMM model
m1 <- readRDS("output/models/MAPE/Fairy Point, Bird Island incubation gamm.rds")

# get smooths
sm <- smooth_estimates(m1$gam, n = 1000) %>%
  add_confint()

# apply exponential to smooths for odds ratios
sm <- sm %>%
  mutate(.estimate = exp(.estimate),
         .lower_ci = exp(.lower_ci),
         .upper_ci = exp(.upper_ci))

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
p7 <- ggplot(cyclones, aes(x = band, fill = as.numeric(OR))) +
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
p7

# plot anticyclones
p8 <- ggplot(anticyclones, aes(x = band, fill = as.numeric(OR))) +
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
p8

# get legend for eddy plots
legend <- get_legend(p7)

# remove legend from eddy plots
p7 <- p7 + theme(legend.position = "none")
p8 <- p8 + theme(legend.position = "none")

# plot eddy plots side by side
odds_ratios <- plot_grid(p7, p8, ncol = 1, rel_heights = c(1, 1)) 
odds_ratios <- plot_grid(odds_ratios, legend, ncol = 2, rel_widths = c(1, 0.25))
odds_ratios <- odds_ratios +
  theme(
    panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent', color = NA), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = NA) #transparent legend panel
  )
odds_ratios

# create eddy shading values for classic GAM plot
ed_vals <- data.frame(type = as.factor(c("anticyclonic core", "anticyclonic periphery", "not eddy",
                                         "cyclonic periphery", "cyclonic core")),
                      x = c(-2.5, -1.5, 0, 1.5, 2.5),
                      xmin = c(-3, -2, -1, 1, 2),
                      xmax = c(-2, -1, 1, 2, 3), 
                      ymin = 0, ymax = 1)

# classic GAM plot
p9 <- ggplot(sm) + 
  geom_tile(data = ed_vals, alpha = 0.4, height = Inf, width = c(1, 1, 2, 1, 1),
            aes(x = x, y = 0.5, fill = type)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = ed2),
              alpha = 0.2) +
  geom_line(aes(x = ed2, y = .estimate), lwd = 0.5) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0,0)) + 
  ylab("Odds Ratio") +
  xlab("Relative Eddy Distance") +
  scale_fill_manual(values = c("red4", "red3",
                                 "steelblue4", "steelblue3",
                                 "white"),
                      breaks = c("anticyclonic core", "anticyclonic periphery",
                                 "cyclonic core", "cyclonic periphery", "not eddy"),
                    name = "Eddy Zone") +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey20")
p9 <- p9 +
  theme(
    panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent', color = NA), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', color = NA) #transparent legend panel
  )

# export
ggsave("output/schematic/classic_gamm.png", p9, bg = "transparent")
ggsave("output/schematic/or_gamm.png", odds_ratios, bg = "transparent")
