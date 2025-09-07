#--------------------------------------------------------
# Plot Eddy Specialisation by CPF Constraints and Species
#--------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

library(tidyverse)
library(cowplot)

# read in data
data <- read.csv("data/tracks/species_site_stage_v3.csv")

# remove NAs in eddy_importance_auger
data <- data %>% filter(eddy_importance_auger != "" &
                          eddy_importance_auger != "remove")

# group by species and cpf 
data <- data %>%
  group_by(species, cpf, eddy_importance_auger) %>%
  summarise(count = n())

# add blank values for each missing combination
species_levels <- unique(data$species)
cpf_levels <- unique(data$cpf)
eddy_levels <- unique(data$eddy_importance_auger)
blank <- expand.grid(species = species_levels, cpf = cpf_levels, eddy_importance_auger = eddy_levels) %>%
  mutate(count = 0)
missing <- blank %>%
  anti_join(data, by = c("species", "cpf", "eddy_importance_auger"))
data <- rbind(data, missing)

# modify factors
data <- data %>%
  mutate(eddy_importance_auger = as.factor(eddy_importance_auger),
         species = as.factor(species)) %>%
  mutate(eddy_importance_auger = fct_relevel(eddy_importance_auger, "specialised", "uncertain", "unspecialised", "not encountered"))

# capitalise
levels(data$eddy_importance_auger) <- c("1. Specialised", "2. Uncertain", "3. Unspecialised", "4. Not Encountered")
levels(data$species) <- c("Ad\u00E9lie", "Chinstrap", "Emperor", "King", "Macaroni")
  
# add numbers for labelling
data <- data %>%
  mutate(number = case_when(
    eddy_importance_auger == "1. Specialised" ~ 1,
    eddy_importance_auger == "2. Uncertain" ~ 2,
    eddy_importance_auger == "3. Unspecialised" ~ 3,
    eddy_importance_auger == "4. Not Encountered" ~ 4
  ))

# if count is 0, set number to NA
data <- data %>% 
  mutate(number = ifelse(count == 0, "", as.character(number)))

# create greyscale legend for eddy_importance_auger
test <- data %>% filter(species == "Ad\u00E9lie")
p1 <- ggplot(test, aes(x = eddy_importance_auger, y = count, fill = eddy_importance_auger)) + geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("grey85", "grey70", "grey55", "grey40"), name = "Eddy Usage") 
p1
legend <- get_legend(p1)

# plot restricted stages
high <- data %>%
  filter(cpf == "high")
p2 <- ggplot(high, aes(x = species, y = count, fill = interaction(species, eddy_importance_auger))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, col = "grey20", lwd = 0.1) +
  geom_text(aes(label=number), position=position_dodge(width=0.8), vjust=-0.25, col = "grey40", size = 3) +
  labs(title = "Highly Constrained Stages", x = "Species", y = "Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#edb191", "#aaf5e8", "#fceaa3", "#e57d7a", "#e6b4ea",
                                "#e69265", "#38e7c9", "#f9dc65", "#da4945", "#d582dc",
                                "#df7339","#118e79", "#f7ce28", "#ba2a25", "#c451ce",
                                "#c65920", "#0a5548", "#d7ae08", "#851e1a", "#76237d"
                                )) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.008, 0.15)) 
plot_grid(p2, plot_grid(NULL, legend, nrow = 2, rel_heights = c(1, 0.8)), rel_widths = c(1, 0.2))
g1 <- plot_grid(p2 + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank()),
                NULL, rel_widths = c(1, 0.2))

# plot medium restricted stages
mid <- data %>%
  filter(cpf == "mid")
p3 <- ggplot(mid, aes(x = species, y = count, fill = interaction(species, eddy_importance_auger))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, col = "grey20", lwd = 0.1) +
  geom_text(aes(label=number), position=position_dodge(width=0.8), vjust=-0.25, col = "grey40", size = 3) +
  labs(title = "Moderately Constrained Stages", x = "Species", y = "Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#edb191", "#aaf5e8", "#fceaa3", "#e57d7a", "#e6b4ea",
                                "#e69265", "#38e7c9", "#f9dc65", "#da4945", "#d582dc",
                                "#df7339","#118e79", "#f7ce28", "#ba2a25", "#c451ce",
                                "#c65920", "#0a5548", "#d7ae08", "#851e1a", "#76237d"
                                )) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.008, 0.15))
g2 <- plot_grid(p3 + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank()),
                plot_grid(NULL, legend, nrow = 2, rel_heights = c(1, 10)), rel_widths = c(1, 0.2))
g2

# plot unrestricted stages
none <- data %>%
  filter(cpf == "none")

p4 <- ggplot(none, aes(x = species, y = count, fill = interaction(species, eddy_importance_auger))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, col = "grey20", lwd = 0.1) +
  geom_text(aes(label=number), position=position_dodge(width=0.8), vjust=-0.25, col = "grey40", size = 3) +
  labs(title = "Unconstrained Stages", x = "Species", y = "Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_fill_manual(values = c("#edb191", "#aaf5e8", "#fceaa3", "#e57d7a", "#e6b4ea",
                                "#e69265", "#38e7c9", "#f9dc65", "#da4945", "#d582dc",
                                "#df7339","#118e79", "#f7ce28", "#ba2a25", "#c451ce",
                                "#c65920", "#0a5548", "#d7ae08", "#851e1a", "#76237d"
                                )) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.008, 0.15))

plot_grid(p4, plot_grid(NULL, legend, nrow = 2, rel_heights = c(1, 0.8)), rel_widths = c(1, 0.2))
g3 <- plot_grid(p4, NULL, rel_widths = c(1, 0.2))

# all stages
g4 <- plot_grid(g1, g2, g3, nrow = 3, rel_heights = c(1, 1, 1.2))
g4

# export
ggsave("output/imagery/bar_importance.png", 
      g4, width = 10, height = 14, dpi = 300)          
