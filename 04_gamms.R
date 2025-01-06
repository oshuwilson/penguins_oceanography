#---------------------------
# ARS Prevalence GAMMs
#---------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(mgcv)
  library(gamm4)
  library(ggside)
  library(patchwork)
}

# define species site and stage
this.species <- "MAPE"
this.site <- "Heard"
this.stage <- "late chick-rearing"

# Load the extracted state tracks
tracks <- readRDS(paste0("output/hmm/hmm_tracks/", this.species, "_", this.site, "_", this.stage, "_tracks_checked.rds"))

# change state to binomial variable
tracks <- tracks %>%
  mutate(bin_state = if_else(state == "ARS", 1, 0))

#resample eddy 0s to values between -0.99 and 0.99 - to avoid misleading trends between -1 and 1
tracks <- tracks %>%
  group_by(row_number()) %>%
  mutate(ed2 = ifelse(eddies == 0, runif(1, -0.99, 0.99), eddies)) %>%
  ungroup()

#create eddie background shading values
ed_vals <- data.frame(type = c("anticyclonic core", "anticyclonic periphery", "not eddy",
                               "cyclonic periphery", "cyclonic core"),
                      x = c(-2.5, -1.5, 0, 1.5, 2.5),
                      xmin = c(-3, -2, -1, 1, 2),
                      xmax = c(-2, -1, 1, 2, 3), 
                      ymin = 0, ymax = 1)

#oceanographic model
m1 <- gamm4(bin_state ~ s(front_freq, k = 5, bs = "ts") + 
              s(eddies, k = 5, bs = "ts") + s(curr, k = 5, bs = "ts"), 
            random = ~(1|individual_id), data = tracks, family = binomial)

# get smooths
sm <- smooth_estimates(m1$gam) %>%
  add_confint()

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

#Front Frequency
p1 <- sm %>%
  filter(.smooth == "s(front_freq)") %>%
  ggplot() + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = front_freq),
              alpha = 0.2
  ) + 
  geom_line(aes(x = front_freq, y = .estimate), lwd = 0.5) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) + 
  ylab("Prevalance of ARS state") +
  xlab("Front Frequency") 
p1 <- p1 + geom_xsidehistogram(data = tracks, aes(x=front_freq), 
                               bins=31, boundary = 0, 
                               fill = "#230493", alpha = 0.2) + 
  theme_ggside_void() + 
  theme(ggside.panel.scale.x = .3) +
  scale_xsidey_continuous(expand = c(0,0))
p1


#Eddies
p2 <- sm %>%
  filter(.smooth == "s(eddies)") %>%
  ggplot() + 
  geom_tile(data = ed_vals, alpha = 0.4, height = Inf, width = c(1, 1, 2, 1, 1),
            aes(x = x, y = 0.5, fill = type)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = eddies),
              alpha = 0.2) + 
  geom_line(aes(x = eddies, y = .estimate), lwd = 0.5) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) + 
  ylab("Prevalance of ARS state") +
  xlab("Relative Eddy Distance")
p2 <- p2 + geom_xsidehistogram(data = tracks, aes(x=ed2), 
                               fill = "#230493", alpha = 0.2,
                               boundary = -3, bins = 60) + 
  theme_ggside_void() + 
  theme(ggside.panel.scale.x = .3,
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical")
p2 <- p2 + scale_fill_manual(values = c("red4", "red3", "white",
                                        "steelblue3", "steelblue4"),
                             breaks = c("anticyclonic core", "anticyclonic periphery",
                                        "not eddy", "cyclonic periphery", "cyclonic core")) +
  scale_xsidey_continuous(expand = c(0,0))
p2

#Current
p3 <- sm %>%
  filter(.smooth == "s(curr)") %>%
  ggplot() + 
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = curr),
              alpha = 0.2) +  
  geom_line(aes(x = curr, y = .estimate), lwd = 0.5) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) + 
  ylab("Prevalance of ARS state") +
  xlab("Current Speed (m/s)")
p3 <- p3 + geom_xsidehistogram(data = tracks, aes(x=curr), 
                               fill = "#230493", alpha = 0.2,
                               boundary = 0, bins = 60) + 
  theme_ggside_void() + 
  theme(ggside.panel.scale.x = .3) + 
  scale_xsidey_continuous(expand = c(0,0))
p3

#plot all together
p1 + p2 + p3 + plot_layout(ncol = 3)
?plot_grid
?plot_layout
# export plots
ggsave(p1, filename = paste0("output/gamms/plots/", this.species, "_", this.site, "_", this.stage, "_front_freq.png"), width = 8, height = 6)
ggsave(p2, filename = paste0("output/gamms/plots/", this.species, "_", this.site, "_", this.stage, "_eddies.png"), width = 8, height = 6)
ggsave(p3, filename = paste0("output/gamms/plots/", this.species, "_", this.site, "_", this.stage, "_current.png"), width = 8, height = 6)

