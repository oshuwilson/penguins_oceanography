##' Generate quality control plots
##'
##' generate a .pdf of SSM fits to RAATD data to aid quality control
##' saves to your existing working directory
##' 
##' @title qcplot
##' 
##' @param fit an object created by using fit_ssm from the aniMotum package
##' @param spp the species' 4-letter abbreviated name
##' @param device_type either GPS, PTT, or GLS
##' 
##' @return a .pdf file of quality control plots
##'
##' @import dplyr
##' @import ggplot2
##' @import aniMotum
##' @import patchwork
##' 
##' @export


qcplot <- function(fit, spp, device_type, step_duration, study_code){
  
  #extract fitted and predicted values
  ggfit <- grab(fit, what = "fitted")
  ggpred <- grab(fit, what = "predicted")
  
  #join with track data to source individual ids
  individuals <- tracks %>% 
    select(id, individual_id, device_id) %>% 
    group_by(id) %>% 
    summarise(individual_id = first(individual_id), device_id = first(device_id))
  ggfit <- ggfit %>% left_join(individuals)
  
  #list of individual/device ids
  inds <- levels(as.factor(ggfit$id))
  
  #setup pdf export
  pdf(file = paste0(spp, "_", device_type, "_", study_code, "qc_plots.pdf"), width = 8, height = 10, pointsize = 16)
  
  #loop plotting over device ids
  for(i in inds){
    #isolate the fitted and prediction for each individual
    this_fit <- ggfit %>% filter(id == i)
    this_pred <- ggpred %>% filter(id == i)
    
    #caculate the standard error for predicted locationsin lat/lon terms
    this_pred$abs_lat <- abs(this_pred$lat)
    40075 * cos(90)/360
    
    #plot 1 - lat/lon for fitted and predicted
    p1 <- ggplot() + geom_point(data=this_fit, aes(x = x, y = y), col = "firebrick", size = 1) +
      geom_point(data=this_pred, aes(x = x, y = y), col = "dodgerblue", size = 1) + 
      geom_path(data=this_pred, aes(x = x, y = y), col = "dodgerblue", lwd = 0.25) +
      theme_bw() + xlab("lat") + ylab("lon") +
      ggtitle(paste("Individual ID: ", this_fit$individual_id, "\n Device ID: ", this_fit$device_id, "\n Device Type: ", device_type, "\n Time Step: ", step_duration, " hour(s)", sep = ""))
    p1 
    
    #plot 2 - time vs lon for fitted and predicted
    p2 <- ggplot() +
      geom_ribbon(data=this_pred, aes(x=date, ymin =x - x.se * 1.96, ymax = x + x.se * 1.96), fill = "dodgerblue4", alpha = 0.3) + 
      geom_point(data=this_fit, aes(x = date, y = x), col = "firebrick", size = 1) +
      geom_rug(data=this_fit, aes(x=date), col = "firebrick", sides = "b") +
      geom_point(data=this_pred, aes(x = date, y = x), col = "dodgerblue", size = 1) + 
      geom_rug(data=this_pred, aes(x=date), col = "dodgerblue", sides = "t") + 
      theme_bw() + ylab("lon")
    p2 
    
    #plot 3 - time vs lat for fitted and predicted
    p3 <- ggplot() +
      geom_ribbon(data=this_pred, aes(x=date, ymin = y - y.se * 1.96, ymax = y + y.se * 1.96), fill = "dodgerblue4", alpha = 0.3) + 
      geom_point(data=this_fit, aes(x = date, y = y), col = "firebrick", size = 1) +
      geom_rug(data=this_fit, aes(x=date), col = "firebrick", sides = "b") +
      geom_point(data=this_pred, aes(x = date, y = y), col = "dodgerblue", size = 1) + 
      geom_rug(data=this_pred, aes(x=date), col = "dodgerblue", sides = "t") + 
      theme_bw() + ylab("lat")
    p3 
    
    #print plot to pdf
    print(p1/p2/p3 + plot_layout(heights = c(5, 2, 2)))
  }
  
  #finish pdf export
  dev.off()
}