#-----------------------------------------------------------
# Plots for incubating Adelie Penguins from Pointe Geologie
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
tracks <- readRDS("output/hmm/hmm_tracks_by_colony/ADPE/Pointe Geologie incubation tracks checked.rds")

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

# only keep dep_years within range of Auger data
deps <- deps %>%
  filter(start_year > 2012 & start_year < 2021)

# get deployment years
dep_years <- deps %>%
  pull(start_year) %>%
  unique() %>%
  sort()


for(z in 1:length(dep_years)){
  
  # for each year
  year <- dep_years[z]
  
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
  eddies <- rast(paste0("E:/Satellite_Data/daily/eddies_auger/eddies_auger_", year, ".nc"))
  eddies <- eddies[[time(eddies) %in% year_dates]]
  
  # crop to extent of tracks
  e <- ext(trax) + c(0.5, 0.5, 0.5, 0.5)
  eddies <- crop(eddies, e)
  
  # remove error layers
  if(year == 2013){
    eddies <- eddies[[time(eddies) != as_date("2013-12-10")]]
  }
  if(year == 2015){
    eddies <- eddies[[time(eddies) != as_date("2015-12-12")]]
  }
  if(year == 2017){
    eddies <- eddies[[!time(eddies) %in% c(as_date("2017-12-13"), as_date("2017-12-22"))]]
  }
  
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
    scale_fill_gradient2(name = "Anticyclone and \nCyclone Prevalence",
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

# select years with biggest sample size (less clutter)
edplots2 <- edplots2[c(3, 5:8)]

# add legend
edplots2[[6]] <- eddy_legend

# plot first 3 years together with legend
edmap <- plot_grid(plotlist = edplots2, nrow = 2)
edmap


#----------------------------------------------------------------------------------
# 2. GAMM circular plots
#----------------------------------------------------------------------------------

# read in GAMM
m1 <- readRDS(paste0("output/models/ADPE/ADPE Pointe Geologie incubation auger gamm.rds"))

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
ggsave("text/draft figs/new/5. Pointe Geologie Incubation.png", grid,
       height = 11, width = 12)


#----------------------------------------------------------------------------------
# 4. Supplementary - Eddies, SIC and Salinity
#----------------------------------------------------------------------------------

# read in depth and crop
depth <- rast("E:/Satellite_Data/static/depth/depth.nc")
depth <- crop(depth, e)

for(z in 1:length(dep_years)){
  
  year <- dep_years[z]
  
  # get individuals belonging to this year
  year_inds <- deps %>% 
    filter(start_year == year) %>%
    pull(individual_id)
  
  # get the track lines for this year
  lines_year <- trax_lines %>%
    filter(individual_id %in% year_inds) 
  
  # limit track to this year
  trax_year <- trax %>%
    filter(individual_id %in% year_inds)
  
  # get dates of these tracks
  year_dates <- as_date(trax_year$date)
  year_dates <- as.POSIXct(year_dates)
  med_date <- as_date(mean(year_dates))
  
  # get SIC for this year
  sic <- rast(paste0("E:/Satellite_Data/daily/sic/sic_", year, ".nc"))
  sic <- sic[[time(sic) == med_date]]

  # crop SIC
  sic <- crop(sic, e)
  
  # if sic is NA, make 0
  sic[is.na(values(sic))] <- 0
  
  # mask out depth NAs
  sic <- mask(sic, depth)
  
  # get a plot of sic for that year
  sicplot <- ggplot() +
    geom_spatraster(data = sic) +
    geom_spatvector(data = crop_coast, fill = "white") +
    geom_spatvector(data = lines_year, col = "grey40", lwd = 0.75) +
    #geom_spatvector(data = trax_year %>% filter(state == "ARS"), aes(col = state), size = 1) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    #scale_color_manual(values = "grey90", name = "", labels = "ARS Events") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle(year) +
    guides(fill = "none", col = "none")
  
  # get sal for this year
  sal <- rast(paste0("E:/Satellite_Data/daily/sal/sal_", year, ".nc"))
  sal <- sal[[time(sal) == med_date]]
  
  # crop sal
  sal <- crop(sal, e)
  
  # get mean sal
  sal <- mean(sal, na.rm = T)
  
  # clamp sal
  sal <- clamp(sal, 33, 35)
  
  # get a plot of sal for that year
  salplot <- ggplot() +
    geom_spatraster(data = sal) +
    geom_spatvector(data = crop_coast, fill = "white") +
    geom_spatvector(data = lines_year, col = "grey40", lwd = 0.75) +
    #geom_spatvector(data = trax_year %>% filter(state == "ARS"), aes(col = state), size = 1) +
    scale_fill_viridis_c(limits = c(33, 35)) +
    #scale_color_manual(values = "grey20", name = "", labels = "ARS Events") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle(year) +
    guides(fill = "none", col = "none")
  
  # # get the eddy plot for the same year
  eddyplot <- edplots[[z]] +
    theme(legend.position = "none")
  
  # plot side by side
  dualplot <- plot_grid(eddyplot, sicplot, salplot, ncol = 3)
  
  # assign to list
  if(z == 1){
    dualplots <- list()
  }
  dualplots[[z]] <- dualplot
}

# limit to example years
dualplots <- dualplots[c(3, 5:8)]

# get legend of sic
siclegplot <- ggplot() +
  geom_spatraster(data = sic * 100) +
  scale_fill_viridis_c(name = "Sea Ice Concentration (%)",
                       limits = c(0, 100))
sic_legend <- get_legend(siclegplot)

# get legend of sal
sallegplot <- ggplot() +
  geom_spatraster(data = sal) +
  scale_fill_viridis_c(name = "Salinity (PSU)",
                       limits = c(33, 35))
sal_legend <- get_legend(sallegplot)

# combine legends
legends <- plot_grid(eddy_legend, sic_legend, sal_legend, ncol = 3)

dualplots[[6]] <- legends

# plot dual plots
dualgrid <- plot_grid(plotlist = dualplots, ncol = 1)
dualgrid

ggsave("text/draft figs/new/S5. Pointe Geologie incubation sic.png",
       dualgrid, width = 15, height = 20)

