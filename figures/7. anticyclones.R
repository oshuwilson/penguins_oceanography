#------------------------------------
# Anticyclones and Nekton
#------------------------------------

# Amanda Bay Fledglings have a very strong tie between anticyclones and all mig_nekton
# Atka Bay Fledglings less so
# Heard chick-rearing Kings not much
# But St Andrews Bay fledglings are very linked

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)
library(terra)
library(tidyterra)
library(mgcv)
library(gratia)

# # read in tracks - either Auger or non
king_tracks <- readRDS("output/hmm/hmm_tracks_by_colony/KIPE/St Andrews Bay, South Georgia late chick-rearing tracks checked.rds")
emp_tracks <- readRDS("output/hmm/hmm_tracks_by_colony/EMPE/Amanda Bay post-breeding tracks checked.rds")

# convert to terra
king_tracks <- king_tracks %>%
  vect(geom = c("x", "y"), crs = "EPSG:6932")
emp_tracks <- emp_tracks %>%
  vect(geom = c("x", "y"), crs = "EPSG:6932")

# read in coastline
coast <- vect("data/coast_vect.RDS")

# read in depth stereographic
depth <- rast("~/OneDrive - University of Southampton/Documents/Predictor Data/depth_stereographic.RDS")

#crop depth to coast extent
depth <- crop(depth, ext(-5790000, 5790000, -5790000, 5790000))

#classify depth into coarser bands
m <- c(0, 1000, 1,
       1000, 2000, 2,
       2000, 3000, 3,
       3000, 4000, 4,
       4000, 5000, 5,
       5000, 6000, 6,
       6000, 7000, 7,
       7000, 8000, 8,
       8000, 9000, 9)
m <- matrix(m, ncol = 3, byrow = TRUE)
depth <- classify(depth, m)

# combine tracks
king_tracks$species <- "King Penguin"
emp_tracks$species <- "Emperor Penguin"
tracks <- rbind(king_tracks, emp_tracks)

# plot all together
map <- ggplot() +
  geom_spatraster(data = depth) +
  geom_spatvector(data = tracks, size = 0.5, aes(col = species)) +
  scale_fill_gradient(na.value = "white",
                      guide = "none",
                      high = "dodgerblue4",
                      low = "lightskyblue1") +
  scale_color_manual(values = c("darkred", "black"), name = "") +
  theme_void() +
  guides(color = guide_legend(override.aes = list(size = 2)))
map

# EMPE tracks only
map_emp <- ggplot() +
  geom_spatraster(data = depth) +
  geom_spatvector(data = emp_tracks, size = 0.5, aes(col = species)) +
  scale_fill_gradient(na.value = "white",
                      guide = "none",
                      high = "dodgerblue4",
                      low = "lightskyblue1") +
  scale_color_manual(values = c("black"), guide = "none") +
  theme_void() +
  xlim(-5e+06, 5e+06) +
  ylim(-5e+06, 5e+06)
map_emp 

# calculate trip lengths
# create lines for each individual
trax <- tracks
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

# project lines to lat/lon
trax_lines <- project(trax_lines, "epsg:4326")

# calculate trip length for each individuals
lengths <- data.frame(track_len = perim(trax_lines), 
                      individual_id = trax_lines$individual_id)
max(lengths$track_len)

# export
#ggsave("output/imagery/fig components/anticyclones_tracks.png", map, width = 10, height = 8, dpi = 300)
ggsave("output/imagery/fig components/anticyclones_tracks_emp.png", map_emp, width = 10, height = 10, dpi = 300)

# 
# # should the updated eddy product be used (depends on timeframe)
# auger <- FALSE
# 
# if(auger == TRUE){
#   tracks$eddies <- tracks$eddies_auger
# }
# 
#
# # convert tracks to terra
# if(auger == TRUE){
#   tracks_terra <- vect(tracks, geom = c("x", "y"), crs = "EPSG:4326")
# } else {
#   tracks_terra <- vect(tracks, geom = c("x", "y"), crs = "EPSG:6932") %>%
#     project("epsg:4326")
# }
# 
# # load extraction functions
# source("code/functions/extraction_functions.R")
# 
# # extract lower_mig_meso_nekton
# tracks_nek <- dynamic_extract("lower_mig_meso_nekton", tracks_terra)
# 
# # plot nekton vs eddies
# ggplot(tracks_nek %>% filter(depth > 0), aes(x=eddies, y = lower_mig_meso_nekton)) + geom_smooth()
# 
# # extract upper_mig_meso_nekton
# tracks_nek <- dynamic_extract("upper_mig_meso_nekton", tracks_nek)
# 
# # plot nekton vs eddies
# ggplot(tracks_nek %>% filter(depth > 0), aes(x=eddies, y = upper_mig_meso_nekton)) + geom_smooth()
# 
# # extract lower_hmig_meso_nekton
# tracks_nek <- dynamic_extract("lower_hmig_meso_nekton", tracks_nek)
# 
# # plot nekton vs eddies
# ggplot(tracks_nek %>% filter(depth > 0), aes(x=eddies, y = lower_hmig_meso_nekton)) + geom_smooth()
# 
# # bank extractions
# saveRDS(tracks_nek, "exploration/stats/anticyclones/SG extractions.rds")

# juvenile emperors feed on lanternfish, which are considered migratory mesopelagic nekton
# let's look at the migratory mesopelagic nekton distribution within different eddies

rm(list=ls())

# read in Amanda Bay and St Andrews Bay fledglings
emps <- readRDS("exploration/stats/anticyclones/Amanda Bay extractions.rds")
kings <- readRDS("exploration/stats/anticyclones/SG extractions.rds")

# # reclass eddy values between -1 and 1 to 0 
# emps <- emps %>%
#   mutate(eddies = ifelse(eddies > -1 & eddies < 1, 0, eddies))
# kings <- kings %>%
#   mutate(eddies = ifelse(eddies > -1 & eddies < 1, 0, eddies))

# # plot eddies vs migratory upper mesopelagic nekton
# emp_up <- ggplot(emps, aes(x=eddies, y = upper_mig_meso_nekton)) +
#   geom_smooth() +
#   labs(x = "Relative Eddy Distance", y = "Migratory Upper Mesopelagic Nekton") +
#   theme_bw()
# 
# king_up <- ggplot(kings, aes(x=eddies, y = upper_mig_meso_nekton)) +
#   geom_smooth() +
#   labs(x = "Relative Eddy Distance", y = "Migratory Upper Mesopelagic Nekton") +
#   theme_bw()
# 
# # plot eddies vs highly migratory lower mesopelagic nekton
# emp_low <- ggplot(emps, aes(x=eddies, y = lower_hmig_meso_nekton)) +
#   geom_smooth() +
#   labs(x = "Relative Eddy Distance", y = "Highly Migratory Lower Mesopelagic Nekton") +
#   theme_bw()
# 
# king_low <- ggplot(kings, aes(x=eddies, y = lower_hmig_meso_nekton)) +
#   geom_smooth() +
#   labs(x = "Relative Eddy Distance", y = "Highly Migratory Lower Mesopelagic Nekton") +
#   theme_bw()
# 
# 
# # modify dataframes to combine upper mig and lower hmig nekton
# emps <- emps %>%
#   as.data.frame() %>%
#   select(eddies, upper_mig_meso_nekton, lower_hmig_meso_nekton) %>%
#   pivot_longer(cols = c(upper_mig_meso_nekton, lower_hmig_meso_nekton), 
#                names_to = "nekton_type", 
#                values_to = "nekton_value")
# kings <- kings %>%
#   as.data.frame() %>%
#   select(eddies, upper_mig_meso_nekton, lower_hmig_meso_nekton) %>%
#   pivot_longer(cols = c(upper_mig_meso_nekton, lower_hmig_meso_nekton), 
#                names_to = "nekton_type", 
#                values_to = "nekton_value")
# 
# # create eddie background shading values
# ed_vals <- data.frame(type = as.factor(c("anticyclonic core", "anticyclonic periphery", "not eddy",
#                                          "cyclonic periphery", "cyclonic core")),
#                       x = c(-2.5, -1.5, 0, 1.5, 2.5),
#                       xmin = c(-3, -2, -1, 1, 2),
#                       xmax = c(-2, -1, 1, 2, 3), 
#                       ymin = 0, ymax = 1)
# 
# # plot upper and lower on same graphs - get linetype legend to display
# emps_both <- ggplot(emps, aes(x=eddies)) +
#   geom_tile(data = ed_vals, alpha = 0.4, height = Inf, width = c(1, 1, 2, 1, 1),
#             aes(x = x, y = 0.5, fill = type)) +
#   geom_smooth(aes(y = nekton_value, linetype = nekton_type), col = "black", fill = "grey40",
#               method = "gam") +
#   labs(x = "Relative Eddy Distance", y = bquote(Nekton ~ Concentration ~ (g ~ m^-2))) +
#   theme_bw() +
#   scale_linetype_manual(values = c("solid", "dashed"), 
#                             labels = c("Highly Migratory Lower Mesopelagic",
#                                    "Upper Migratory Mesopelagic"),
#                         name = "Type of Nekton") +
#   ggtitle("Fledgling Emperor Penguins from Amanda Bay, Enderby Land") +
#   scale_y_continuous(expand = c(0,0)) +
#   coord_cartesian(ylim = c(0, 4)) +
#   scale_fill_manual(values = c("red4", "red3",
#                                "steelblue4", "steelblue3",
#                                "grey80"),
#                     breaks = c("anticyclonic core", "anticyclonic periphery",
#                                "cyclonic core", "cyclonic periphery", "not eddy"),
#                     labels = c("Anticyclonic Core", "Anticyclonic Periphery",
#                                "Cyclonic Core", "Cyclonic Periphery", "Not Eddy"),
#                     name = "Eddy Zone") +
#   scale_x_continuous(expand = c(0,0), limits = c(-3, 3)) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# emps_both
# 
# king_both <- ggplot(kings, aes(x=eddies)) +
#   geom_tile(data = ed_vals, alpha = 0.4, height = Inf, width = c(1, 1, 2, 1, 1),
#             aes(x = x, y = 0.5, fill = type)) +
#   geom_smooth(aes(y = nekton_value, linetype = nekton_type), col = "black", fill = "grey40",
#               method = "gam") +
#   labs(x = "Relative Eddy Distance", y = bquote(Nekton ~ Concentration ~ (g ~ m^-2))) +
#   theme_bw() +
#   scale_linetype_manual(values = c("solid", "dashed"), 
#                             labels = c("Highly Migratory \nLower Mesopelagic Nekton",
#                                    "Upper Migratory \nMesopelagic Nekton"),
#                         name = "Type of Nekton") +
#   ggtitle("Fledgling King Penguins from St Andrews Bay, South Georgia") +
#   scale_y_continuous(expand = c(0,0)) +
#   coord_cartesian(ylim = c(1, 1.6)) +
#   scale_fill_manual(values = c("red4", "red3",
#                                "steelblue4", "steelblue3",
#                                "grey80"),
#                     breaks = c("anticyclonic core", "anticyclonic periphery",
#                                "cyclonic core", "cyclonic periphery", "not eddy"),
#                     labels = c("Anticyclonic Core", "Anticyclonic Periphery",
#                                "Cyclonic Core", "Cyclonic Periphery", "Not Eddy"),
#                     name = "Eddy Zone") +
#   scale_x_continuous(expand = c(0,0), limits = c(-3, 3)) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# king_both
# 
# # plot together
# library(cowplot)
# 
# # get legend from emps_both
# legend <- get_legend(emps_both)
# 
# # remove legend from both plots
# emps_both <- emps_both + theme(legend.position = "none")
# king_both <- king_both + theme(legend.position = "none")
# 
# # remove xlab from emps_both
# emps_both <- emps_both + labs(x = NULL)
# 
# # combine plots
# combined_plot <- plot_grid(emps_both, king_both, ncol = 1, align = "hv")
# combined_plot
# 
# # add legend to combined plot
# combined_plot_with_legend <- plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(1, 0.4))
# combined_plot_with_legend
# 
# # export the combined plot
# ggsave("text/draft figs/7. anticyclones.png", 
#        plot = combined_plot_with_legend, 
#        width = 9, height = 10, dpi = 300)

# reclass eddy values between -1 and 1 to 0
emps <- emps %>%
   mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))
# kings <- kings %>%
#   mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))
  
# fit GAMs to extractions
m1 <- gam(upper_mig_meso_nekton ~ s(ed2), data = emps)
plot(m1)
m2 <- gam(lower_hmig_meso_nekton ~ s(ed2), data = emps)
plot(m2)
m3 <- gam(upper_mig_meso_nekton ~ s(ed2), data = kings)
plot(m3)
m4 <- gam(lower_hmig_meso_nekton ~ s(ed2), data = kings)
plot(m4)

# get smooths
sm1 <- smooth_estimates(m1, n = 1000) %>%
  add_confint() %>%
  add_constant(coef(m1)[1])
sm2 <- smooth_estimates(m2, n = 1000) %>%
  add_confint() %>%
  add_constant(coef(m1)[1])
sm3 <- smooth_estimates(m3, n = 1000) %>%
  add_confint()  
sm4 <- smooth_estimates(m4, n = 1000) %>%
  add_confint()

# average odds ratios over .5 intervals
sm1 <- sm1 %>%
  mutate(ed2 = cut(ed2, breaks = c(-3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.125, 
                                   0.125, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3))) %>%
  group_by(ed2) %>%
  summarise(.estimate = mean(.estimate),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci))
sm2 <- sm2 %>%
  mutate(ed2 = cut(ed2, breaks = c(-3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.125, 
                                   0.125, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3))) %>%
  group_by(ed2) %>%
  summarise(.estimate = mean(.estimate),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci))

sm3 <- sm3 %>%
  mutate(ed2 = cut(ed2, breaks = c(-3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.125, 
                                   0.125, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3))) %>%
  group_by(ed2) %>%
  summarise(.estimate = mean(.estimate),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci))

sm4 <- sm4 %>%
  mutate(ed2 = cut(ed2, breaks = c(-3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.125, 
                                   0.125, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3))) %>%
  group_by(ed2) %>%
  summarise(.estimate = mean(.estimate),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci))

# get upper and lower limits from cuts
sm1 <- sm1 %>%
  mutate(x_tmp = str_sub(ed2, 2, -2)) %>%
  separate(x_tmp, c("min", "max"), sep = ",") %>%
  mutate_at(c("min", "max"), as.numeric)
sm2 <- sm2 %>%
  mutate(x_tmp = str_sub(ed2, 2, -2)) %>%
  separate(x_tmp, c("min", "max"), sep = ",") %>%
  mutate_at(c("min", "max"), as.numeric)
sm3 <- sm3 %>%
  mutate(x_tmp = str_sub(ed2, 2, -2)) %>%
  separate(x_tmp, c("min", "max"), sep = ",") %>%
  mutate_at(c("min", "max"), as.numeric)
sm4 <- sm4 %>%
  mutate(x_tmp = str_sub(ed2, 2, -2)) %>%
  separate(x_tmp, c("min", "max"), sep = ",") %>%
  mutate_at(c("min", "max"), as.numeric)

# remove two categories with no associated eddy region
sm1 <- sm1 %>% 
  filter(!ed2 %in% c("(-1,-0.125]", "(0.125,1]"))
sm2 <- sm2 %>% 
  filter(!ed2 %in% c("(-1,-0.125]", "(0.125,1]"))
sm3 <- sm3 %>%
  filter(!ed2 %in% c("(-1,-0.125]", "(0.125,1]"))
sm4 <- sm4 %>%
  filter(!ed2 %in% c("(-1,-0.125]", "(0.125,1]"))

# rename odds ratio
sm1 <- sm1 %>%
  rename(OR = .estimate)
sm2 <- sm2 %>%
  rename(OR = .estimate)
sm3 <- sm3 %>%
  rename(OR = .estimate)
sm4 <- sm4 %>%
  rename(OR = .estimate)

# get estimate for no eddy zone
no_eddy1 <- sm1 %>%
  filter(ed2 == "(-0.125,0.125]") %>%
  select(OR) %>%
  pull()
no_eddy2 <- sm2 %>%
  filter(ed2 == "(-0.125,0.125]") %>%
  select(OR) %>%
  pull()
no_eddy3 <- sm3 %>%
  filter(ed2 == "(-0.125,0.125]") %>%
  select(OR) %>%
  pull()
no_eddy4 <- sm4 %>%
  filter(ed2 == "(-0.125,0.125]") %>%
  select(OR) %>%
  pull()

# subtract no_eddy stat from all odds ratios
sm1 <- sm1 %>%
  mutate(OR = OR - no_eddy1)
sm2 <- sm2 %>%
  mutate(OR = OR - no_eddy2)
sm3 <- sm3 %>%
  mutate(OR = OR - no_eddy3)
sm4 <- sm4 %>%
  mutate(OR = OR - no_eddy4)

# create cyclonic dataframe
cyclones1 <- sm1 %>% 
  filter(min >= 1) %>%
  arrange(desc(min))
cyclones2 <- sm2 %>% 
  filter(min >= 1) %>%
  arrange(desc(min))
cyclones3 <- sm3 %>%
  filter(min >= 1) %>%
  arrange(desc(min))
cyclones4 <- sm4 %>%
  filter(min >= 1) %>%
  arrange(desc(min))

# create factor for plotting
cyclones1 <- cyclones1 %>%
  mutate(band = factor(1:nrow(cyclones1)))
cyclones2 <- cyclones2 %>%
  mutate(band = factor(1:nrow(cyclones2)))
cyclones3 <- cyclones3 %>%
  mutate(band = factor(1:nrow(cyclones3)))
cyclones4 <- cyclones4 %>%
  mutate(band = factor(1:nrow(cyclones4)))

# create anticyclonic dataframe
anticyclones1 <- sm1 %>% 
  filter(min <= -1) %>%
  arrange(min)
anticyclones2 <- sm2 %>%
  filter(min <= -1) %>%
  arrange(min)
anticyclones3 <- sm3 %>%
  filter(min <= -1) %>%
  arrange(min)
anticyclones4 <- sm4 %>%
  filter(min <= -1) %>%
  arrange(min)

# create factor for plotting
anticyclones1 <- anticyclones1 %>%
  mutate(band = factor(1:nrow(anticyclones1)))
anticyclones2 <- anticyclones2 %>%
  mutate(band = factor(1:nrow(anticyclones2)))
anticyclones3 <- anticyclones3 %>%
  mutate(band = factor(1:nrow(anticyclones3)))
anticyclones4 <- anticyclones4 %>%
  mutate(band = factor(1:nrow(anticyclones4)))

# combine emp data frames
emp_cyclones <- cyclones1 %>%
  mutate(OR = OR + cyclones2$OR)
emp_anticyclones <- anticyclones1 %>%
  mutate(OR = OR + anticyclones2$OR)
king_cyclones <- cyclones3 %>%
  mutate(OR = OR + cyclones4$OR)
king_anticyclones <- anticyclones3 %>%
  mutate(OR = OR + anticyclones4$OR)

# get min and max values
emp_min <- min(c(emp_cyclones$OR, emp_anticyclones$OR), na.rm = TRUE)
emp_max <- max(c(emp_cyclones$OR, emp_anticyclones$OR), na.rm = TRUE)
king_min <- min(c(king_cyclones$OR, king_anticyclones$OR), na.rm = TRUE)
king_max <- max(c(king_cyclones$OR, king_anticyclones$OR), na.rm = TRUE)
# min_val1 <- min(sm1$OR, na.rm = TRUE)
# max_val1 <- max(sm1$OR, na.rm = TRUE)
# min_val2 <- min(sm2$OR, na.rm = TRUE)
# max_val2 <- max(sm2$OR, na.rm = TRUE)
# min_val3 <- min(sm3$OR, na.rm = TRUE)
# max_val3 <- max(sm3$OR, na.rm = TRUE)
# min_val4 <- min(sm4$OR, na.rm = TRUE)
# max_val4 <- max(sm4$OR, na.rm = TRUE)

# take biggest of both absolute values
emp_abs_val <- max(abs(emp_min), abs(emp_max))
king_abs_val <- max(abs(king_min), abs(king_max))
# abs_val1 <- max(abs(min_val1), abs(max_val1))
# abs_val2 <- max(abs(min_val2), abs(max_val2))
# abs_val3 <- max(abs(min_val3), abs(max_val3))
# abs_val4 <- max(abs(min_val4), abs(max_val4))

# change min and max to absolute
emp_min <- -emp_abs_val
emp_max <- emp_abs_val
king_min <- -king_abs_val
king_max <- king_abs_val
# min_val1 <- -abs_val1
# max_val1 <- abs_val1
# min_val2 <- -abs_val2
# max_val2 <- abs_val2
# min_val3 <- -abs_val3
# max_val3 <- abs_val3
# min_val4 <- -abs_val4
# max_val4 <- abs_val4

# delta
d <- "Delta"

emp1 <- ggplot(emp_cyclones, aes(x = band, fill = as.numeric(OR))) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 8.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 0,
                       limits = c(emp_min, emp_max),
                       na.value = "grey") +
  ggtitle("Cyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5)) +
  #labs(fill = expression(Delta~" Nekton Concentration (g m "*-2*")")) +
  labs(fill = bquote(.(sym(d)) ~ Nekton ~ Concentration ~ (g ~ m^-2))) +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")) +
  guides(fill = guide_colorbar(title.position = "top"))

emp2 <- ggplot(emp_anticyclones, aes(x = band, fill = as.numeric(OR))) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 8.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 0,
                       limits = c(emp_min, emp_max),
                       na.value = "grey") +
  ggtitle("Anticyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5)) +
  #labs(fill = expression(Delta~" Nekton Concentration (g m "*-2*")")) +
  labs(fill = bquote(.(sym(d)) ~ Nekton ~ Concentration ~ (g ~ m^-2))) +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")) +
  guides(fill = guide_colorbar(title.position = "top"))

ggplot(king_cyclones, aes(x = band, fill = as.numeric(OR))) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 8.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 0,
                       limits = c(king_min, king_max),
                       na.value = "grey") +
  ggtitle("Cyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5)) +
  #labs(fill = expression(Delta~" Nekton Concentration (g m "*-2*")")) +
  labs(fill = bquote(.(sym(d)) ~ Nekton ~ Concentration ~ (g ~ m^-2))) +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")) +
  guides(fill = guide_colorbar(title.position = "top"))

ggplot(king_anticyclones, aes(x = band, fill = as.numeric(OR))) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 8.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 0,
                       limits = c(king_min, king_max),
                       na.value = "grey") +
  ggtitle("Anticyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5)) +
  #labs(fill = expression(Delta~" Nekton Concentration (g m "*-2*")")) +
  labs(fill = bquote(.(sym(d)) ~ Nekton ~ Concentration ~ (g ~ m^-2))) +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")) +
  guides(fill = guide_colorbar(title.position = "top"))


# check Emperor plots
emp1 
emp2

# plot together
library(cowplot)

eddies <- plot_grid(emp1, emp2)


# export
ggsave("output/imagery/fig components/anticyclones_emp_cyclones.png", emp1, width = 8, height = 6, dpi = 300)
ggsave("output/imagery/fig components/anticyclones_emp_anticyclones.png", emp2, width = 8, height = 6, dpi = 300)
ggsave("text/draft figs/new/S6. Nekton Concentrations.png", width = 12, height = 8)
