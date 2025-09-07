#-----------------------------------------------
# Make colony summary plots for each case study
#-----------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 02")

{
  library(tidyverse)
  library(terra)
  library(tidyterra)
  library(suncalc)
  library(gratia)
  library(gamm4)
  library(cowplot)
}

# read in species-site-stage metadata
srs <- read.csv("data/tracks/species_site_stage_v2.csv")
# FIX i = 31, 33, 34
i <- 53

for(i in 35:nrow(srs)){
  
  # define species, site, and stage
  this.species <- srs$species[i]
  this.site <- srs$site[i]
  this.stage <- srs$stage[i]
  area <- srs$island[i]
  auger <- srs$auger[i]
  
  # coastline
  coast <- rnaturalearth::ne_countries(scale = 10, returnclass = "sv")
  
  
  #-----------------------------------------------------------------------------
  # 1. Eddy maps by deployment year
  #-----------------------------------------------------------------------------
  
  # read in tracks
  tracks <- readRDS(paste0("output/hmm/hmm_tracks_by_colony/", this.species, "/", this.site, " ", this.stage, " tracks checked.rds"))
  
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
  
  # get deployment years
  deps <- tracks %>%
    group_by(individual_id) %>%
    summarise(start_date = min(as_date(date))) %>%
    mutate(start_year = year(start_date))
  dep_years <- deps %>%
    pull(start_year) %>%
    unique() %>%
    sort()
  
  # remove any outside satellite data range
  dep_years <- dep_years[dep_years < 2022]
  if(auger == "yes"){
    dep_years <- dep_years[dep_years > 2012]
  } else {
    dep_years <- dep_years[dep_years > 1992]
  }
  
  for(z in 1:length(dep_years)){
    
    # for each year
    year <- dep_years[z]
    yearn <- year + 1
    
    # get tracks belonging to this year
    year_inds <- deps %>% 
      filter(start_year == year) %>%
      pull(individual_id)
    tracks_year <- trax %>%
      filter(individual_id %in% year_inds)
    
    # get dates of these tracks
    year_dates <- unique(as_date(tracks_year$date))
    
    # if fewer than 3 days, skip
    if(length(year_dates) < 3){
      next
    }
    
    # read in eddies for this year
    if(auger == "yes"){
      eddies <- rast(paste0("E:/Satellite_Data/daily/eddies_auger/eddies_auger_", year, ".nc"))
    } else {
      eddies <- rast(paste0("E:/Satellite_Data/daily/eddies/eddies_", year, ".nc"))
    }
    eddies <- eddies[[time(eddies) %in% year_dates]]
    
    # if any dates fall into the following year, also read in the following year
    if(yearn %in% year(year_dates)){
      if(auger == "yes"){
        eddiesn <- rast(paste0("E:/Satellite_Data/daily/eddies_auger/eddies_auger_", yearn, ".nc"))
      } else {
        eddiesn <- rast(paste0("E:/Satellite_Data/daily/eddies/eddies_", yearn, ".nc"))
      }
      eddiesn <- eddiesn[[time(eddiesn) %in% year_dates]]
      eddies <- c(eddies, eddiesn)
    }
    
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
    crop_coast <- crop(coast, ext(trax) + c(0.5, 0.5, 0.5, 0.5))
    
    # plot this year's tracks
    p1 <- ggplot() +
      geom_spatraster(data = P) +
      geom_spatvector(data = crop_coast, fill = "white") +
      geom_spatvector(data = lines_year, col = "grey40", lwd = 0.75) +
      geom_spatvector(data = tracks_year %>% filter(state == "ARS"), aes(col = state), size = 1) +
      scale_fill_gradient2(name = "Anticyclone and \ncyclone prevalence",
                           low = "#F08080", high = "steelblue4", limits = c(-1, 1)) +
      scale_color_manual(values = "grey20", name = "", labels = "ARS Events") +
      theme_bw() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      ggtitle(year)
    print(p1)
    
    # store plot in list
    if(!exists("edplots")){
      edplots <- list()
    }
    edplots[[z]] <- p1
    
  }
  
  
  #----------------------------------------------------------------------------------
  # 2. GAMM circular plots
  #----------------------------------------------------------------------------------
  
  # load in GAMM model
  if(auger == "yes"){
    m1 <- readRDS(paste0("output/models/", this.species, "/", this.species, " ", this.site, " ", this.stage, " auger gamm.rds"))
  } else { 
    try(m1 <- readRDS(paste0("output/models/", this.species, "/", this.site, " ", this.stage, " gamm.rds")))
  }
  
  # if model exists (i.e. eddies encountered)
  if(exists("m1")){
    
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
    
  } else { #if no eddies detected in tracks
    
    # read in blank dataframes
    cyclones <- readRDS("output/imagery/colony odds ratios/cyclones blank.rds")
    anticyclones <- readRDS("output/imagery/colony odds ratios/anticyclones blank.rds")
    bg <- data.frame(ed2 = "(-0.125,0.125]", OR = as.numeric(NA), 
                     .lower_ci = as.numeric(NA), .upper_ci = as.numeric(NA),
                     x = 1, y = 1)
    
  }
  
  
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
  
  
  #-------------------------------------------------------------------------------
  # 3. GAMM statistics
  #-------------------------------------------------------------------------------
  
  # read in GAMM stats
  stats <- readRDS("output/models/model_stats.rds")
  
  # filter to this case study
  stats <- stats %>%
    filter(species == this.species &
             site == this.site &
             stage == this.stage)
  
  # rename variables
  stats <- stats %>%
    mutate(variable = case_when(
      variable == "ed2" ~ "Eddies",
      variable == "depth" ~ "Depth",
      variable == "curr" ~ "Current",
      variable == "sic" ~ "Sea Ice"
    ))
  
  # mutate values of 0 to < 2e-16
  stats <- stats %>%
    mutate(`p-value` = ifelse(`p-value` == 0, "< 2e-16", paste0("= ", `p-value`)))
  
  # format label text
  label_lines <- apply(stats, 1, function(row) {
    sprintf("%s: χ² = %.1f, p %s", row["variable"], as.numeric(row["Chi.sq"]), (row["p-value"]))
  })
  label_text <- paste(label_lines, collapse = "\n")
  label_text
  
  # create label
  plot_label <- ggdraw() +
    draw_label(label_text)
  
  
  #-------------------------------------------------------------------------------
  # 4. Bring it all together
  
  # get legend from the first eddy plot
  eddy_legend <- get_legend(edplots[[1]])
  
  # remove legend from all eddy plots
  edplots2 <- lapply(edplots, function(x){
    x <- x + theme(legend.position = "none")
  })
  
  # if more than 3 plots, keep three at random
  if(length(edplots2) > 3){
    edplots3 <- sample(edplots2, 3)
  } else {
    edplots3 <- edplots2
  }
  
  # plot first 3 years together with legend
  edmap <- plot_grid(plotlist = edplots3, nrow = 1)
  edmap <- plot_grid(edmap, eddy_legend, ncol = 2, rel_widths = c(1, 0.2))
  
  # stack with odds ratios
  final_plot <- plot_grid(edmap, odds_ratios, plot_label, ncol = 1, rel_heights = c(1, 0.6, 0.3))
  print(final_plot)
  
  # export
  ggsave(paste0("output/imagery/colony summary plots new/", this.species, "/", this.site, " ", this.stage, " overview.png"), 
         final_plot, height = 12, width = 13, units = "in")
  
  # remove conditional items
  rm(m1, edplots)
}

