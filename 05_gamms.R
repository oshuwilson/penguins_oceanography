#---------------------------
# ARS Prevalence GAMMs
#---------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(mgcv)
  library(gamm4)
  library(gratia)
}

# define species site and stage
this.species <- "KIPE"
this.site <- "Macquarie"
this.stage <- "incubation"

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

p1 <- draw(m1$gam, select = "s(front_freq)", rug = F, constant = coef(m1$gam)[1],
           fun = inv_link(m1$gam), ci_col = "black")

p2 <- draw(m1$gam, select = "s(eddies)", rug = F, constant = coef(m1$gam)[1],
           fun = inv_link(m1$gam), ci_col = "black")

p3 <- draw(m1$gam, select = "s(curr)", rug = F, constant = coef(m1$gam)[1],
           fun = inv_link(m1$gam), ci_col = "black")


p1 <- p1 + theme_bw() +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  ggtitle(paste0("Prevalance of hidden states for ", this.species, " at ", this.site, " (", this.stage, ")")) +
  ylab("Prevalance of ARS state") +
  xlab("Front Frequency")
p1

p2 <- p2 + theme_bw() +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  ggtitle(paste0("Prevalance of hidden states for ", this.species, " at ", this.site, " (", this.stage, ")")) +
  ylab("Prevalance of ARS state") +
  xlab("Eddies") + 
  geom_tile(data = ed_vals, alpha = 0.4, height = Inf, width = c(1, 1, 2, 1, 1),
            aes(x = x, y = 0.5, fill = type)) + 
  scale_fill_manual(values = c("red4", "red3", "steelblue4", "steelblue3", "white")) + 
  theme(legend.position = "bottom")
p2

p3 <- p3 + theme_bw() +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  ggtitle(paste0("Prevalance of hidden states for ", this.species, " at ", this.site, " (", this.stage, ")")) +
  ylab("Prevalance of ARS state") +
  xlab("Current Speed (m/s)")
p3

# export plots
ggsave(p1, filename = paste0("output/gamms/plots/", this.species, "_", this.site, "_", this.stage, "_front_freq.png"), width = 8, height = 6)
ggsave(p2, filename = paste0("output/gamms/plots/", this.species, "_", this.site, "_", this.stage, "_eddies.png"), width = 8, height = 6)
ggsave(p3, filename = paste0("output/gamms/plots/", this.species, "_", this.site, "_", this.stage, "_current.png"), width = 8, height = 6)

