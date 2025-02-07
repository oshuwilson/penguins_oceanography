#---------------------------------------------------
# Assess eddy usage by feeding time or no. of events
#---------------------------------------------------

#remove nighttime points (and large error?) before GAMM
#add depth (and curr/SIC?) to GAMMs
#if nothing available from a category in background, change to greyed out

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(gamm4)
  library(gratia)
}

# example population 
this.species <- "KIPE"
this.site <- "Hound Bay, South Georgia"
this.stage <- "chick-rearing"

# load in hmm checked tracks
tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))


# 1. Process Data

# filter to ARS only
ars <- tracks %>% filter(state == "ARS")

# resample non-eddies to between -1 and 1 to avoid artificial inflation
ars <- ars %>%
  mutate(ed2 = ifelse(eddies > -1 & eddies < 1, runif(n(), -0.99999, 0.99999), eddies))
ars <- ars %>%
  mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))

# histogram of foraging
ggplot(ars, aes(x = ed2)) + geom_histogram(boundary = -3, breaks = seq(-3, 3, 0.1))

# read in background samples
back <- readRDS(paste0("exploration/background_OR/", this.species, "/", this.site, " ", this.stage, " background.rds"))

# only keep ARS-corresponding samples
back <- back %>% 
  as.data.frame(geom = "XY")

# resample non-eddies in background
back <- back %>%
  mutate(ed2 = ifelse(eddies > -1 & eddies < 1, runif(n(), -0.99999, 0.99999), eddies))
back <- back %>%
  mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))

# create binary presence/absence cols
ars$pa <- 1
back$pa <- 0


# 2. GAMM

# join datasets together
data <- ars %>% 
  select(individual_id, ed2, depth, curr, pa) %>%
  rbind(select(back, individual_id, ed2, depth, curr, pa))

# run GAMM on eddies
m1 <- gamm4(pa ~ s(ed2, bs = "ts"),
            random = ~(1|individual_id), family = binomial, data = data)

# get smooths
sm <- smooth_estimates(m1$gam) %>%
  add_confint()
sm

# add constant
constant <- coef(m1$gam)[1]
sm <- sm %>%
  mutate(.estimate = .estimate + constant,
         .lower_ci = .lower_ci + constant,
         .upper_ci = .upper_ci + constant)

# transform the smooths
backtrans <- inv_link(m1$gam)
sm <- sm %>%
  mutate(.estimate = backtrans(.estimate),
         .lower_ci = backtrans(.lower_ci),
         .upper_ci = backtrans(.upper_ci))

# plot
ggplot(sm) + geom_line(aes(x = ed2, y = .estimate)) + geom_ribbon(aes(x = ed2, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) + theme_minimal() + labs(x = "Eddy usage", y = "Presence probability")


# 3. Odds Ratios

# reset smooths
sm <- smooth_estimates(m4$gam, n = 1000) %>%
  add_confint()

# apply exponential to smooths for odds ratios
sm <- sm %>%
  mutate(.estimate = exp(.estimate),
         .lower_ci = exp(.lower_ci),
         .upper_ci = exp(.upper_ci))

# plot
ggplot(sm) + 
  geom_line(aes(x = ed2, y = .estimate)) + 
  geom_ribbon(aes(x = ed2, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) + 
  theme_minimal() + 
  labs(x = "Relative Eddy Distance", y = "Odds Ratio") + 
  scale_y_continuous(limits = c(0, 3))
ggsave("exploration/background_OR/gamm_tests/eddy_depth_curr_gamm.png")

# average odds ratios over .5 intervals
sm2 <- sm %>%
  mutate(ed2 = cut(ed2, breaks = c(-3, -2.5, -2, -1.5, -1, -0.25, 
                                   0.25, 1, 1.5, 2, 2.5, 3))) %>%
  group_by(ed2) %>%
  summarise(.estimate = mean(.estimate),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci))

# remove two categories with no associated eddy region
sm2 <- sm2[-c(5,7),]

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

# plot cyclones
ggplot(cyclones, aes(x = band, fill = OR)) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 2.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 1,
                       name = "Odds Ratio", limits = c(0, 3),
                       labels = c(0, 1, 2, "3+")) +
  ggtitle("Cyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))
ggsave("exploration/background_OR/gamm_tests/eddy_depth_curr_cyclones.png")

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

# plot anticyclones
ggplot(anticyclones, aes(x = band, fill = OR)) +
  geom_bar(width = 1) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 2.5, color = "black", lwd = 1) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", lwd = 0.65) +
  geom_vline(xintercept = 4.5, color = "black", lwd = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 1,
                       name = "Odds Ratio", limits = c(0, 3),
                       labels = c(0, 1, 2, "3+")) +
  ggtitle("Anticyclones") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

ggsave("exploration/background_OR/gamm_tests/eddy_depth_curr_anticyclones.png")
