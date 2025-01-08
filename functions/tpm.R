##' tpm
##'
##' Plot the transition probability matrix for a given covariate 
##' for each individual and the mean trend line
##'
##' @title tpm
##' 
##' @param model hmm model produced with hmmTMB::hmm$fit()
##' @param var the covariate to be plotted 
##' @param inds a list of each individual ID code 
##'
##' @import hmmTMB
##' @import foreach
##' @import doParallel
##' @import tidyverse
##' 
##' @export

tpm <- function(model, var, inds){
  
  #set number of cores for parallelisation
  ncores <- length(inds)
  if(ncores > 10){
    ncores <- 10
  }
  
  #enable parallelisation 
  cl <- makePSOCKcluster(10)
  registerDoParallel(cl)
  
  #foreach loop over each individual
  alldata <- foreach(i = inds, .combine = rbind, .packages = c("hmmTMB", "ggplot2", "dplyr")) %dopar% {
    
    #set individual to this individual
    covs <- list(individual_id = i)
    
    #plot covariates
    pn <- model$plot("tpm", var = var, i = 2, j = 1, covs = covs)
    
    #extract geom_line data
    plotdata <- pn$data %>%
      select(var, prob) %>% #keep x and y columns
      mutate(ind = i) #add individual ID
    
    #combine this data
    return(plotdata)
    
  }
  
  #stop parallelisation
  stopCluster(cl)
  
  #calculate mean trend
  mean <- alldata %>%
    group_by(var) %>%
    summarise(prob = mean(prob, na.rm = T))
  
  #plot for all individuals
  p1 <- ggplot(alldata, aes(x = var, y = prob)) +
    geom_line(aes(group = ind), linewidth = 0.7, color = "lightgrey") + 
    geom_line(data = mean, aes(x = var, y = prob), color = "darkred", linewidth = 1) +
    ylab("Transition Probability from Transit to ARS") +
    xlab(var) + 
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, max(alldata$prob, na.rm = T)), expand = c(0, 0)) +
    ggtitle(" ") 
  
  return(p1)
  
}