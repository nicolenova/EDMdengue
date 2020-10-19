######################################
# Functions for Forecasting Analyses #
#                                    #
# nicole.nova@stanford.edu           #
######################################


# Forecasts to investigate the predictability of the model and its attractor 
forecast.train <- function(df, ts.dim, forecast.horizon){
  
  # Length of new dataframe
  dfl <- length(df$date)
  
  # Edge/outer part of a simplex
  n_start <- 2
  n_end <- 2
  
  # Creating simplexes using data from half a year, a year, and 1.5 years ago
  previous.simplexes <- c(seq(forecast.horizon, 20), seq(50, 59), seq(74, 88))
  
  # Initialize forecast results for each forecasting week
  ds.f <- list()
  
  for (nn in n_start:n_end) { 
    
    for (sf in previous.simplexes) { # Varying over simplex forecasts
      
      ssr.cases <- block_lnlp(df, columns = ts.dim, 
                              target_column = "cases", stats_only = F, 
                              method = "simplex", exclusion_radius = NULL, tp =  sf, 
                              num_neighbors = nn) 
        
      stats <- ssr.cases[[1]]$stats
      ssr.cases <- ssr.cases[[1]]$model_output
      ssr.cases.obs <- ssr.cases$obs*sd(PR.orig$cases) + mean(PR.orig$cases)
      ssr.cases.pred <- ssr.cases$pred*sd(PR.orig$cases) + mean(PR.orig$cases) # Predicted time series
      
      preds <- data.frame(ssr.cases.obs, ssr.cases.pred)
      
      ds.f[[which(previous.simplexes == sf)]] <- preds  
    }
  }
  
  ts.length <- dfl - (dfl - length(18:900))
  
  ds.wk <- list()
  wk <- 0
  for (i in (ts.length-1):0) {
    wk <- wk+1
    ds.time <- NULL
    for (sf in previous.simplexes) {
      ds.time <- c(ds.time, ds.f[[which(previous.simplexes == sf)]]$ssr.cases.pred[dfl-sf-i])
    }
    ds.wk[[wk]] <- ds.time
  }
  
  # Now each week has a forecast distribution so we fit to skewed normal
  low.q <- NULL
  up.q <- NULL
  mean.sn <- NULL
  
  time <- 1:ts.length
  
  preds <- ds.wk
  
  for (t in time) {
    # Lower quantile
    low.q <- c(low.q, qsnorm(0.05, 
                             mean = snormFit(preds[[t]])$par[[1]], 
                             sd = snormFit(preds[[t]])$par[[2]], 
                             if (snormFit(preds[[t]])$par[[3]] > 5) {
                               xi = 5  # Keep large numbers but avoid infinity
                             } else {
                               xi = snormFit(preds[[t]])$par[[3]]
                             }) )
    
    # Upper quantile
    up.q <- c(up.q, qsnorm(0.95, 
                           mean = snormFit(preds[[t]])$par[[1]], 
                           sd = snormFit(preds[[t]])$par[[2]], 
                           if (snormFit(preds[[t]])$par[[3]] > 5) {
                             xi = 5  # Keep large numbers but avoid infinity
                           } else {
                             xi = snormFit(preds[[t]])$par[[3]]
                           }) )
    
    # Mean
    mean.sn <- c(mean.sn, snormFit(preds[[t]])$par[[1]] )
  }
  
  start.wk <- dfl - ts.length + 1
  
  low <- data.frame(x = start.wk:dfl, y = low.q)
  low$y[low$y < 0] <- 0
  #low <- as.data.frame(supsmu(low$x, low$y, bass = 1))   # Optional smoothing

  up <- data.frame(x = start.wk:dfl, y = up.q)
  up$y[up$y < 0] <- 0
  #up <- as.data.frame(supsmu(up$x, up$y, bass = 1))      # Optional smoothing

  m <- data.frame(x = start.wk:dfl, y = mean.sn)
  m$y[m$y < 0] <- 0
  
  # Combine to one dataframe
  forecasts <- data.frame(week = start.wk:dfl, 
                          date = PR$date[start.wk:dfl],
                          obs = ssr.cases.obs[18:900],
                          low.q = low$y,
                          mean.pred = m$y,
                          up.q = up$y)
  
  return(forecasts)
}






# Out of sample forecasts to test model performance
forecast.test <- function(df, ts.dim, forecast.horizon, previous.simplexes){
  
  quarter <- 13 # weeks
  train.length <- 988 # length of the original training dataset in weeks
  
  # Create new dataframes with added quarters sequentially throughout the testing seasons
  for (i in 1:16) {
    end <- train.length + quarter*i
    assign(paste0("PR.", i), PR[1:end,])
  }
  
  # Edge/outer part of a simplex
  n_start <- 2
  n_end <- 2
  
  # Initialize time series for concatenating out of sample forecasts
  forecasts <- data.frame()
  
  # # Specify which dataset to use for forecasting and update the attractor for new forecasts
  # for (i in 1:16) {
  # 
  #   df <- eval(parse(text = paste0("PR.", i)))
  # # There is no difference in forecast performance when adding new data, weekly, quarterly or yearly
  # # but there is a big difference in computing power
    
  for (dfNumber in seq(4, 16, 4)) {

    df <- eval(parse(text = paste0("PR.", dfNumber)))
    
    # Length of new dataframe
    dfl <- length(df$date)
    
    # Initialize forecast results for each forecasting week
    ds.f <- list()
  
    for (sf in previous.simplexes) { # Varying over simplex forecasts
      counter <- 0
      for (nn in n_start:n_end) {
        counter <- counter+1
        
        ssr.cases <- block_lnlp(df, columns = ts.dim, 
                                target_column = "cases", stats_only = F, 
                                method = "simplex", exclusion_radius = NULL, tp =  sf, 
                                num_neighbors = nn) 
        
        stats <- ssr.cases[[1]]$stats
        ssr.cases <- ssr.cases[[1]]$model_output
        ssr.cases.obs <- ssr.cases$obs*sd(PR.orig$cases) + mean(PR.orig$cases)
        ssr.cases.pred <- ssr.cases$pred*sd(PR.orig$cases) + mean(PR.orig$cases) # Predicted time series
        
        if (counter == 1) {
          preds <- data.frame(obs = ssr.cases.obs, pred = ssr.cases.pred)
        } else {
          preds[ , ncol(preds) + 1] <- ssr.cases.pred            # Append new column
          colnames(preds)[ncol(preds)] <- paste0("pred.nn_", nn)
        }
      }
      
      ds.f[[which(previous.simplexes == sf)]] <- preds  
    }
    
    # For testing time series we only need data starting from week 989
    ts.length <- dfl - train.length
    
    ds.wk <- list()
    wk <- 0
    for (i in (ts.length-1):0) {
      wk <- wk+1
      ds.time <- NULL
      for (sf in previous.simplexes) {
        ndf <- ds.f[[which(previous.simplexes == sf)]]
        ds.time <- c(ds.time, as.numeric(ndf[dfl-sf-i, 2:ncol(ndf)]))
      }
      ds.wk[[wk]] <- ds.time
    }
    
    # Now each week has a forecast distribution so we fit to skewed normal
    low.q <- NULL
    up.q <- NULL
    mean.sn <- NULL
    
    time <- 1:ts.length
    
    preds <- ds.wk
    
    for (t in time) {
      # Lower quantile
      low.q <- c(low.q, qsnorm(0.05, 
                               mean = snormFit(preds[[t]])$par[[1]], 
                               sd = snormFit(preds[[t]])$par[[2]], 
                               if (snormFit(preds[[t]])$par[[3]] > 5) {
                                 xi = 5   # Keep large numbers but avoid infinity
                               } else {
                                 xi = snormFit(preds[[t]])$par[[3]]
                               }) )
      
      # Upper quantile
      up.q <- c(up.q, qsnorm(0.95, 
                             mean = snormFit(preds[[t]])$par[[1]], 
                             sd = snormFit(preds[[t]])$par[[2]], 
                             if (snormFit(preds[[t]])$par[[3]] > 5) {
                               xi = 5   # Keep large numbers but avoid infinity
                             } else {
                               xi = snormFit(preds[[t]])$par[[3]]
                             }) )
      
      # Mean
      mean.sn <- c(mean.sn, snormFit(preds[[t]])$par[[1]] )
    }
    
    start.wk <- dfl - ts.length + 1
    
    low <- data.frame(x = start.wk:dfl, y = low.q)
    low$y[low$y < 0] <- 0
    #low <- as.data.frame(supsmu(low$x, low$y, bass = 1))   # Optional smoothing
    
    up <- data.frame(x = start.wk:dfl, y = up.q)
    up$y[up$y < 0] <- 0
    #up <- as.data.frame(supsmu(up$x, up$y, bass = 1))      # Optional smoothing

    m <- data.frame(x = start.wk:dfl, y = mean.sn)
    m$y[m$y < 0] <- 0
    
    # Combine to one dataframe
    last.sf <- previous.simplexes[length(previous.simplexes)]
    forecasts.segment <- data.frame(week = start.wk:dfl, 
                                    date = df$date[start.wk:dfl],
                                    obs = ssr.cases.obs[(start.wk-last.sf):(dfl-last.sf)],
                                    low.q = low$y,
                                    mean.pred = m$y,
                                    up.q = up$y)
    
    #df.start <- nrow(forecasts.segment)-quarter+1
    df.start <- nrow(forecasts.segment)-quarter*4+1
    
    # Add each segment to the full forecast for all testing years
    forecasts <- rbind(forecasts, forecasts.segment[df.start:nrow(forecasts.segment), ] )
  }
  
  return(forecasts)
}







# Plotting the forecasts
gen.forecast.plot <- function(model, mainTitle, timelinePlacement){
  
  # Obtain correlation coefficient (prediction skill) 
  c <- cor.test(model$obs, model$mean.pred)
  rhoval <- round(c$estimate, digits = 4)
  pval <- round(c$p.value, digits = 10)
  
  df <- model
  
  # Set colors
  col_obs <- "black" 
  col_pred <- "#00BFC4" 
  col_bound <- rgb(col2rgb(col_pred)[1]/255, col2rgb(col_pred)[2]/255, col2rgb(col_pred)[3]/255, 0.4)
  
  ggplot(df, aes(x = date)) +
    geom_line(aes(y = low.q), color = "white", size = 1) +
    geom_line(aes(y = up.q), color = "white", size = 1) +
    geom_ribbon(aes(ymin = low.q, ymax = up.q), fill = col_bound) +
    geom_line(aes(y = mean.pred), color = col_pred, size = 1.1) +
    geom_line(aes(y = obs), color = col_obs) +
    xlab("Year") +
    ylab("Incidence") +
    ylim(0, 461) +
    annotate("text", size = 6, x = as.Date(timelinePlacement), y = 400, label = bquote(rho == .(rhoval) )) +
    #annotate("text", size = 6, x = as.Date(timelinePlacement), y = 320, label = bquote("p < 0.001")) +
    ggtitle(mainTitle) +
    theme(text = element_text(size = 25), 
          axis.text = element_text(size = 18),
          plot.title = element_text(size = 18))
     
}







# Save subplots
plot.forecast <- function(model, subplot_name, mainTitle, trainORtest, timelinePlacement){
  
  subplot <- gen.forecast.plot(model, mainTitle, timelinePlacement)
  
  print(subplot)
  file_name = paste("../output/fig_4/", trainORtest, "/", subplot_name, ".pdf", sep="")
  dev.copy(pdf, file=file_name, width = 8, height = 3)
  dev.off()
}







# END OF SCRIPT