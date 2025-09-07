#-----------------------------------------------------------
# Plots for fledgling Emperor Penguins from Amanda Bay
#-----------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02/")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(suncalc)
  library(gratia)
  library(gamm4)
  library(cowplot)
}

# coastline
coast <- rnaturalearth::ne_countries(scale = 10, returnclass = "sv")

# read in tracks
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/EMPE/Amanda Bay post-breeding tracks checked.rds")

# convert to terra
trax <- tracks %>%
  vect(geom = c("x", "y"), crs = "epsg:6932") %>%
  project("epsg:4326")

# reconvert to dataframe
tracks <- trax %>% as.data.frame(geom = "XY") %>%
  rename(lon = x, lat = y, datetime = date) %>%
  mutate(date = as.Date(datetime))

#get dawn times
tracks$dawn <- getSunlightTimes(data = tracks,
                                keep = c("dawn"), tz = "UTC") %>%
  pull(dawn)

#get dusk times
tracks$dusk <- getSunlightTimes(data = tracks,
                                keep = c("dusk"), tz = "UTC") %>%
  pull(dusk)

#if point is between dawn and dusk, then label as daytime
tracks <- tracks %>%
  mutate(timeofday = ifelse(datetime >= dawn & datetime <= dusk, "day", "night"))

# if an ARS event happens during the night, reclass as other behaviour
tracks <- tracks %>%
  mutate(state = case_when(
    state == "ARS" & timeofday == "night" ~ "other",
    is.na(timeofday) ~ state,
    TRUE ~ state))

# recreate trax
trax <- tracks %>%
  vect(geom = c("lon", "lat"), crs = "epsg:4326")

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

# list dates
dates <- unique(as_date(tracks$date))

# list years
years <- unique(year(dates)) %>% sort()

# get deployment info
deps <- tracks %>%
  group_by(individual_id) %>%
  summarise(start_date = min(as_date(date))) %>%
  mutate(start_year = year(start_date))

# get deployment years
dep_years <- deps %>%
  pull(start_year) %>%
  unique() %>%
  sort()

for(z in 1:length(dep_years)){
  
  # for each year
  year <- dep_years[z]
  yearn <- year + 1
  
  # get individuals belonging to this year
  year_inds <- deps %>% 
    filter(start_year == year) %>%
    pull(individual_id)
  
  # get tracks belonging to this year
  tracks_year <- trax %>%
    filter(individual_id %in% year_inds)
  
  # get dates of these tracks
  year_dates <- unique(as_date(tracks_year$date))
  
  # if fewer than 3 days, skip
  if(length(year_dates) < 3){
    next
  }
  
  # read in eddies for these dates
  eddies <- rast(paste0("E:/Satellite_Data/daily/eddies/eddies_", year, ".nc"))
  eddies <- eddies[[time(eddies) %in% year_dates]]
  eddiesn <- rast(paste0("E:/Satellite_Data/daily/eddies/eddies_", yearn, ".nc"))
  eddiesn <- eddiesn[[time(eddiesn) %in% year_dates]]
  eddies <- c(eddies, eddiesn)
  
  # crop to extent of tracks
  e <- ext(trax) + c(0.5, 0.5, 0.5, 0.5)
  eddies <- crop(eddies, e)
  
  # classify eddies as anticyclone (below -1), cyclone (above 1) or non-eddy
  mat1 <- matrix(c(-3, -1, -1,
                   -1, 1, 0,
                   1, 3, 1),
                 nrow = 3,
                 byrow = T)
  
  eds <- classify(eddies, mat1, include.lowest = TRUE, right = TRUE)
  
  # number of days
  N <- nlyr(eds)
  
  # calculate number of anticyclone and cyclone days
  if(N > 1){
    C <- app(eds, fun = function(x) sum(x == 1, na.rm = TRUE)/N)
    A <- app(eds, fun = function(x) sum(x == -1, na.rm = TRUE)/N)
  } else {
    next
  }
  
  # total eddy days ratio
  R <- sum(C, A)
  
  # eddy presence polarity
  P <- (C-A)
  P[is.nan(P)] <- 0
  
  # get the track lines for this year
  lines_year <- trax_lines %>%
    filter(individual_id %in% year_inds) 
  
  # crop coast to tracks
  crop_coast <- crop(coast, e)
  
  # plot this year's tracks
  p1 <- ggplot() +
    geom_spatraster(data = P) +
    geom_spatvector(data = crop_coast, fill = "white") +
    geom_spatvector(data = lines_year, col = "grey40", lwd = 0.75) +
    geom_spatvector(data = tracks_year %>% filter(state == "ARS"), aes(col = state), size = 0.5) +
    scale_fill_gradient2(name = "Anticyclone and \ncyclone prevalence",
                         low = "#F08080", high = "steelblue4", limits = c(-1, 1)) +
    scale_color_manual(values = "grey20", name = "", labels = "ARS Events") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle(year)
  print(p1)
  
  # store plot in list
  if(z == 1){
    edplots <- list()
  }
  edplots[[z]] <- p1
}


# get legend from the first eddy plot
eddy_legend <- get_legend(edplots[[1]])

# remove legend from all eddy plots
edplots2 <- lapply(edplots, function(x){
  x <- x + theme(legend.position = "none")
})

# plot first 3 years together with legend
edmap <- plot_grid(plotlist = edplots2, nrow = 2)

# plot with the legend
edmap <- plot_grid(edmap, eddy_legend, rel_widths = c(0.8, 0.2))
edmap


#----------------------------------------------------------------------------------
# 2. GAMM circular plots
#----------------------------------------------------------------------------------

# read in GAMM
m1 <- readRDS(paste0("output/models/EMPE/Amanda Bay post-breeding gamm.rds"))

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

# get odds ratio for background values
bg <- sm2 %>%
  filter(ed2 == "(-0.125,0.125]") %>%
  mutate(x = 1, y = 1)


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

# plot background values
p4 <- ggplot(bg, aes(x = x, y = y, fill = OR)) +
  geom_tile() +
  coord_fixed() + 
  scale_fill_gradient2(low = "steelblue4", high = "darkred", midpoint = 1,
                       name = "Odds Ratio", limits = c(0, 3),
                       labels = c(0, 1, 2, "3+"),
                       na.value = "grey") +
  theme_void() +  
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black")) +
  ggtitle("Outside\nEddies") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5))

# get legend for eddy plots
legend <- get_legend(p3)

# stack legend and p4
mid <- plot_grid(legend, p4, ncol = 1)

# remove legend from eddy plots
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

# plot eddy plots side by side
odds_ratios <- plot_grid(p3, mid, p2, ncol = 3, rel_widths = c(1, 0.2, 1)) 
odds_ratios


#----------------------------------------------------------------------------------
# 3. Plot all together
#----------------------------------------------------------------------------------

# plot eddies and odds ratios
grid <- plot_grid(edmap, odds_ratios, ncol = 1, rel_heights = c(0.7, 0.3),
                  labels = "auto")
grid

# export
ggsave("text/draft figs/new/6. Amanda Bay Fledglings.png", grid, 
       height = 12, width = 11)


#----------------------------------------------------------------------------------
# 4. Micronekton
#----------------------------------------------------------------------------------

# nekton extracted data
emps <- readRDS("exploration/stats/anticyclones/Amanda Bay extractions.rds")

# reclass eddy values between -1 and 1 to 0
emps <- emps %>%
  mutate(ed2 = ifelse(eddies > -1 & eddies < 1, 0, eddies))

# combine micronekton
emps <- emps %>%
  mutate(micronek = upper_mig_meso_nekton + 
           lower_hmig_meso_nekton)

# fit GAM to extractions
m1 <- gam(micronek ~ s(ed2), data = emps)
plot(m1)

# get smooths
sm1 <- smooth_estimates(m1, n = 1000) %>%
  add_confint() %>%
  add_constant(coef(m1)[1])

# plot smooths
ggplot(sm1, aes(x = ed2, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
              fill = "grey") +
  geom_line() +
  theme_bw()

max(sm1$.estimate)
sm1 %>%
  filter(ed2 < 0.01 & ed2 > 0) %>%
  pull(.estimate)
sm1 %>%
  filter(ed2 > 2.8) %>%
  pull(.estimate)
summary(m1)
