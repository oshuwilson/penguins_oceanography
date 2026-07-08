#-------------------------------------------------------------------------------
# Test whether eddy amplitude/age is greater in those used by penguins
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

library(tidyverse)
library(lme4)
library(lmerTest)
library(cowplot)
library(emmeans)

# list all the in-depth case studies
eddy_cases <- list.files("output/eddy attributes/", full.names = TRUE, pattern = "*.rds")

# choose case study (1 to 5)
i <- 1
print(eddy_cases[i])

# read in the eddy attributes
eddy_attributes <- readRDS(eddy_cases[i])

# change radius to km and amplitude to cm
eddy_attributes <- eddy_attributes %>%
  mutate(eddy_radius = eddy_radius / 1000,
         eddy_amplitude = eddy_amplitude * 100)

# change pa to presence and background
eddy_attributes <- eddy_attributes %>%
  mutate(pa = ifelse(pa == 1, "ARS", "Background"))

# capitalise eddy type
eddy_attributes <- eddy_attributes %>%
  mutate(eddy_type = case_when(
    eddy_type == "anticyclone" ~ "Anticyclones",
    eddy_type == "cyclone" ~ "Cyclones",
    TRUE ~ eddy_type
  ))

# combine pa with eddy_type
eddy_attributes <- eddy_attributes %>%
  mutate(group = paste(eddy_type, pa, sep = " "))

# plot eddy age
p1 <- ggplot(eddy_attributes, aes(x = eddy_maturity)) +
  geom_density(aes(fill = group), alpha = 0.5, color = NA) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Density") +
  xlab("Eddy Maturity") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5),
        strip.text = element_text(size = 12)) +
  scale_fill_manual(name = "Location Type", 
                    values = c("Anticyclones ARS" = "red3", 
                               "Anticyclones Background" = "grey30",
                               "Cyclones ARS" = "steelblue1",
                               "Cyclones Background" = "grey30")) +
  facet_wrap(~eddy_type, ncol = 1, scales = "free_y")
p1

# get average age of eddies used by penguins vs background
eddy_attributes %>%
  group_by(group) %>%
  summarise(median_age = median(eddy_maturity, na.rm = TRUE),
            mean_age = mean(eddy_maturity, na.rm = TRUE),
            sd_age = sd(eddy_maturity, na.rm = TRUE))

# mixed effects model to compare whether eddy age is greater in penguin eddies
m1 <- glmer(eddy_maturity ~ group + (group | individual_id) + (1 | individual_id),
            data = eddy_attributes,
            family = binomial)
summary(m1)
anova(m1)
em1 <- emmeans(m1, pairwise ~ group, type = "response", adjust = "none")
summary(em1)


# plot eddy amplitude
p2 <- ggplot(eddy_attributes, aes(x = eddy_amplitude)) +
  geom_density(aes(fill = group), alpha = 0.5, color = NA) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Density") +
  xlab("Eddy Amplitude (cm)") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5),
        strip.text = element_text(size = 12)) +
  scale_fill_manual(name = "Location Type", 
                    values = c("Anticyclones ARS" = "red3", 
                               "Anticyclones Background" = "grey30",
                               "Cyclones ARS" = "steelblue1",
                               "Cyclones Background" = "grey30")) +
  facet_wrap(~eddy_type, ncol = 1, scales = "free_y")
p2

# get average amplitude of eddies used by penguins vs background
eddy_attributes %>%
  group_by(group) %>%
  summarise(median_amp = median(eddy_amplitude, na.rm = TRUE),
            mean_amp = mean(eddy_amplitude, na.rm = TRUE),
            sd_amp = sd(eddy_amplitude, na.rm = TRUE))


# mixed effects model to compare whether eddy amplitude is greater in penguin eddies
m2 <- glmer(eddy_amplitude ~ group + (group | individual_id) + (1 | individual_id),
            data = eddy_attributes %>% filter(eddy_amplitude != 0),
            family = Gamma(link = "log"))
summary(m2)
em2 <- emmeans(m2, pairwise ~ group, type = "response", adjust = "none")
summary(em2)

# calculate eddy intensity
eddy_attributes <- eddy_attributes %>%
  mutate(intensity = eddy_amplitude / eddy_radius)

# plot eddy intensity
p4 <- ggplot(eddy_attributes, aes(x = intensity)) +
  geom_density(aes(fill = group), alpha = 0.5, color = NA) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Density") +
  xlab("Eddy Intensity (cm/km)") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5),
        strip.text = element_text(size = 12)) +
  scale_fill_manual(name = "Location Type", 
                    values = c("Anticyclones ARS" = "red3", 
                               "Anticyclones Background" = "grey30",
                               "Cyclones ARS" = "steelblue1",
                               "Cyclones Background" = "grey30")) +
  facet_wrap(~eddy_type, ncol = 1, scales = "free_y")
p4

# get average intensity of eddies used by penguins vs background
eddy_attributes %>%
  group_by(group) %>%
  summarise(median_intensity = median(intensity, na.rm = TRUE),
            mean_intensity = mean(intensity, na.rm = TRUE),
            sd_intensity = sd(intensity, na.rm = TRUE))

# mixed effects model to compare whether eddy intensity is greater in penguin eddies
m_intensity <- glmer(intensity ~ group + (group | individual_id) + (1 | individual_id),
            data = eddy_attributes %>% filter(intensity != 0),
            family = Gamma(link = "log"))
summary(m_intensity)
em_intensity <- emmeans(m_intensity, pairwise ~ group, adjust = "none")
summary(em_intensity)

# # model eke, age and amplitude together
# eddy_attributes$pa <- ifelse(eddy_attributes$pa == "ARS", 1, 0)
# 
# m3 <- glmer(
#   pa ~ log(eddy_age) + log(eddy_amplitude) + log(eke) + (1 | individual_id),
#   family = binomial,
#   data = eddy_attributes %>% filter(eddy_amplitude != 0 & eddy_age != 0)
# )
# summary(m3)

# get legend
legend <- get_legend(p2)

# plot all plots together
grid <- plot_grid(p1 + theme(legend.position = "none"),
          p2 + theme(legend.position = "none") + ylab(""),
          p4 + theme(legend.position = "none") + ylab(""),
          legend, ncol = 4,
          rel_widths = c(1, 1, 1, 0.5),
          align = "v", axis = "tb")
grid + ggview::canvas(width = 12, height = 6)

# export grid
ggsave(paste0("output/eddy attributes/plots/eddy_attributes_case_", i, ".png"), grid, width = 12, height = 6)
saveRDS(grid, paste0("output/eddy attributes/plots/eddy_attributes_case_", i, ".rds"))
