################################################
# Functions for Forecasting Challenge Analyses #
#                                              #
# nicole.nova@stanford.edu                     #
################################################


# Out of sample forecasts
# Generating forecasts while updating the model/attractor using only specified past training data
forecast.test <- function(PR, train.length, ts.dim, target){
  
  # Edge/outer part of a simplex
  n_start <- 2
  n_end <- 2
  
  # Creating simplexes using data from half a year, a year, and 1.5 years ago
  previous.simplexes <- c(seq(1, 20), seq(50, 59), seq(74, 88))
  
  # Generate multiple datasets with varying lengths
  quarter <- 13 # weeks
  train.length <- train.length # length of the original training dataset in weeks
  
  # Create new dataframes with added quarters sequentially thoughout the testing seasons
  for (i in 0:16) {
    end <- train.length + quarter*i
    assign(paste0("PR.", i), PR[1:end, ])
  }
  
  # Set parameters for length of data on the first week in each season (1-4)
  s1 <- length(PR.0$cases) + 1
  s2 <- length(PR.4$cases) + 1
  s3 <- length(PR.8$cases) + 1
  s4 <- length(PR.12$cases) + 1
  
  # Specify which dataset to use for forecasting and update the attractor for new forecasts
  for (dfNumber in seq(4, 16, 4)) { 
    
    df <- eval(parse(text = paste0("PR.", dfNumber))) 
    
    # Length of new dataframe
    dfl <- length(df$date)
    
    # Season number
    season <- dfNumber/4
    
    # Obtain the denormalized observations for a particular season
    obs.s <- df$cases*sd(PR.orig$cases) + mean(PR.orig$cases)
    start.season <- eval(parse(text = paste0("s", season)))
    obs.s <- obs.s[start.season:length(df$cases)]
    
    # Initialize forecast results for each forecasting week f_0, f_4, ..., f_24
    sample1 <- NULL
    sample2 <- NULL
    sample3 <- NULL
    sample4 <- NULL
    sample5 <- NULL
    sample6 <- NULL
    sample7 <- NULL
    
    for (nn in n_start:n_end) { 
      
      for (sf in previous.simplexes) { # Varying over simplex forecasts starting from closest week to season
      
        # Initialize the forecasts for each seasons on weeks 0, 4, ..., 24
        f_0 <- NULL
        f_4 <- NULL
        f_8 <- NULL
        f_12 <- NULL
        f_16 <- NULL
        f_20 <- NULL
        f_24 <- NULL
        
        # Add new observations for f_4, f_8, f_12, f_16, f_20, f_24
        f_4 <- c(f_4, obs.s[1:4])
        f_8 <- c(f_8, obs.s[1:8])
        f_12 <- c(f_12, obs.s[1:12])
        f_16 <- c(f_16, obs.s[1:16])
        f_20 <- c(f_20, obs.s[1:20])
        f_24 <- c(f_24, obs.s[1:24])
        
        s <- start.season - sf
        
        for (forecast.horizon in sf:(51+sf)) {   # For each forecast horizon: tp = 1, 2, .., 52 
          
          ssr.cases <- block_lnlp(df, columns = ts.dim.cases.drivers, 
                                  target_column = "cases", stats_only = F, 
                                  method = "simplex", exclusion_radius = NULL, tp =  forecast.horizon, 
                                  num_neighbors = nn) 
          
          stats <- ssr.cases[[1]]$stats
          ssr.cases <- ssr.cases[[1]]$model_output
          ssr.cases.obs <- ssr.cases$obs*sd(PR.orig$cases) + mean(PR.orig$cases)
          ssr.cases.pred <- ssr.cases$pred*sd(PR.orig$cases) + mean(PR.orig$cases)
          
          f_0 <- c(f_0, ssr.cases.pred[s])
          f_4 <- c(f_4, ssr.cases.pred[s+4])
          f_8 <- c(f_8, ssr.cases.pred[s+8])
          f_12 <- c(f_12, ssr.cases.pred[s+12]) 
          f_16 <- c(f_16, ssr.cases.pred[s+16])
          f_20 <- c(f_20, ssr.cases.pred[s+20]) 
          f_24 <- c(f_24, ssr.cases.pred[s+24]) 
        }
        
        e <- 52 # weeks in a year
        f <- list(f_0[1:e], f_4[1:e], f_8[1:e], f_12[1:e], f_16[1:e], f_20[1:e], f_24[1:e]) 
        
        if (target == 0) {
          # These are the n forecasts for each season
          sample1 <- append(sample1, list(f[[1]]))
          sample2 <- append(sample2, list(f[[2]]))
          sample3 <- append(sample3, list(f[[3]]))
          sample4 <- append(sample4, list(f[[4]]))
          sample5 <- append(sample5, list(f[[5]]))
          sample6 <- append(sample6, list(f[[6]]))
          sample7 <- append(sample7, list(f[[7]]))
        } else if (target == 1) {
          # These are the estimates of peak week for each n forecast
          sample1 <- c(sample1, which.max(f[[1]]))
          sample2 <- c(sample2, which.max(f[[2]]))
          sample3 <- c(sample3, which.max(f[[3]]))
          sample4 <- c(sample4, which.max(f[[4]]))
          sample5 <- c(sample5, which.max(f[[5]]))
          sample6 <- c(sample6, which.max(f[[6]]))
          sample7 <- c(sample7, which.max(f[[7]]))
        } else if (target == 2) {
          # These are the estimates of peak incidence for each n forecast
          sample1 <- c(sample1, max(f[[1]]))
          sample2 <- c(sample2, max(f[[2]]))
          sample3 <- c(sample3, max(f[[3]]))
          sample4 <- c(sample4, max(f[[4]]))
          sample5 <- c(sample5, max(f[[5]]))
          sample6 <- c(sample6, max(f[[6]]))
          sample7 <- c(sample7, max(f[[7]]))
        } else {
          # These are the estimates of seasonal incidence for each n forecast
          sample1 <- c(sample1, sum(f[[1]]))
          sample2 <- c(sample2, sum(f[[2]]))
          sample3 <- c(sample3, sum(f[[3]]))
          sample4 <- c(sample4, sum(f[[4]]))
          sample5 <- c(sample5, sum(f[[5]]))
          sample6 <- c(sample6, sum(f[[6]]))
          sample7 <- c(sample7, sum(f[[7]]))
        }

      }
    }
    
    # Collect forecasts for each forecasting week for each season (1-4)
    assign(paste0("sample.", season), list(sample1, sample2, sample3, sample4, sample5, sample6, sample7) ) 
    
  } # End of for loop for each season
  
  forecast.seasons <- list(sample.1, sample.2, sample.3, sample.4)
  
  return(forecast.seasons)
}







# Calculating the peak score from distribution
calc.score.peak.wk <- function(forecast.seasons, seasonNumber, obs.sx) {
  
  sfMax <- 45 # length(previous.simplexes)
  
  # Week where the true peak is
  obs <- which.max(obs.sx) 
  
  sx.p_i <- NULL
  p <- NULL
  
  for (i in 1:7) { # For each forecasting week 0, 4, ..., 24
    # Generate the distribution based on varying forecasts
    preds <- forecast.seasons[[seasonNumber]][[i]][1:sfMax]
    
    if (var(preds) == 0) {
      p <- c(p, 1)
    } else {
      ds <- NULL
      for (wk in 1:52) {
        ds <- c(ds, dsnorm(wk, mean = snormFit(preds)$par[[1]], 
                           sd = snormFit(preds)$par[[2]],
                           xi = snormFit(preds)$par[[3]]) )
      }
      
      # The distribution from all sf forecasts
      plot(dsnorm(seq(0, 52, by = 1), mean = snormFit(preds)$par[[1]],
                  sd = snormFit(preds)$par[[2]],
                  xi = snormFit(preds)$par[[3]]), type="l", ylab = "")
    }
    
    p <- ds
    p_i <- p[obs]
    log(p_i)
    sx.p_i <- c(sx.p_i, log(p_i))
  }
  
  score.sx <- mean(sx.p_i)
  
  return(score.sx)
}








# Calculating the score for peak incidence or seasonal incidence
calc.score.peak.tot.inc <- function(forecast.seasons, seasonNumber, opt.obs.sx, binMin, binMax, maxInc) {
  
  sfMax <- 45 # length(previous.simplexes)
  
  p <- NULL
  
  for (i in 1:7) {  # For each forecasting week 0, 4, ..., 24
    
    preds <- forecast.seasons[[seasonNumber]][[i]][1:sfMax]
    
    if (var(preds) == 0) {
      p <- c(p, 1)
    } else {
      p <- c(p, psnorm(binMax, mean = snormFit(preds)$par[[1]], 
                       sd = snormFit(preds)$par[[2]], 
                       xi = snormFit(preds)$par[[3]]) 
             - psnorm(binMin, mean = snormFit(preds)$par[[1]], 
                      sd = snormFit(preds)$par[[2]], 
                      xi = snormFit(preds)$par[[3]])) 
      
      plot(dsnorm(seq(0, maxInc, by = 1), mean = snormFit(preds)$par[[1]],
                  sd = snormFit(preds)$par[[2]],
                  xi = snormFit(preds)$par[[3]]), type="l", ylab = "")
      abline(v = c(binMin, opt.obs.sx, binMax))
    }
  }
  
  inc.score <- mean(log(p))
  
  return(inc.score)
}







# Plotting forecasting figures
plot.forecasts <- function(forecast.seasons, seasonNumber, obs.sx, maxInc, fig) {
  
  sfMax <- 45 # length(previous.simplexes)
  time <- 1:52
  
  for (i in 1:7) {  # For each forecasting week 0, 4, ..., 24
    
    wk <- seq(0, 24, 4)[[i]]
  
    preds <- list()
    for (t in time) {
      time.ds <- NULL  # Distribution across sf for a time point
      for (sf in 1:sfMax) { # For each time point obtain ds metrics across sf
        time.ds <- c(time.ds, forecast.seasons[[seasonNumber]][[i]][[sf]][[t]])
      }
      preds[[t]] <- time.ds
    }
    
    # Values of a target
    low.q <- NULL
    mean.sn <- NULL
    med.q <- NULL
    up.q <- NULL
    
    for (t in time) {
      # Lower quantile
      low.q <- c(low.q, qsnorm(0.025, 
                               mean = snormFit(preds[[t]])$par[[1]], 
                               sd = snormFit(preds[[t]])$par[[2]], 
                               if (snormFit(preds[[t]])$par[[3]] > 5) {
                                 xi = 5   # Keep large numbers but avoid infinity
                               } else {
                                 xi = snormFit(preds[[t]])$par[[3]]
                               }) )
      
      # Mean
      mean.sn <- c(mean.sn, snormFit(preds[[t]])$par[[1]] )
      
      # Upper quantile
      up.q <- c(up.q, qsnorm(0.975, 
                             mean = snormFit(preds[[t]])$par[[1]], 
                             sd = snormFit(preds[[t]])$par[[2]], 
                             xi = snormFit(preds[[t]])$par[[3]]) )
      
    }
    
    
    if (is.na(low.q[1]) == T) {
      low.q[1] <- 0
    } else {
      low.q[1] <- low.q[1]
    }
    
    if (is.na(low.q[length(low.q)]) == T) {
      low.q[length(low.q)] <- 0
    } else {
      low.q[length(low.q)] <- low.q[length(low.q)]
    }
    
    
    if (is.na(up.q[1]) == T) {
      up.q[1] <- 0
    } else {
      up.q[1] <- up.q[1]
    }
    
    if (is.na(up.q[length(up.q)]) == T) {
      up.q[length(up.q)] <- 0
    } else {
      up.q[length(up.q)] <- up.q[length(up.q)]
    }
    
    # Add observations for a particular forecasting week
    obs <- data.frame(x = time[1:wk], y = obs.sx[1:wk])
      
    low <- data.frame(x = time[wk:52], y = low.q[wk:52])
    slow <- as.data.frame(supsmu(low$x, low$y, bass = 2))
    slow$y[slow$y < 0] <- 0
    slow <- slow[-1,]
    lower <- rbind(obs, slow)
    
    up <- data.frame(x = time[wk:52], y = up.q[wk:52])
    sup <- as.data.frame(supsmu(up$x, up$y, bass = 2))
    sup$y[sup$y < 0] <- 0
    sup <- sup[-1,]
    upper <- rbind(obs, sup)

    # Average forecast
    m <- data.frame(x = time[wk:52], y = mean.sn[wk:52])
    sm <- as.data.frame(supsmu(m$x, m$y, bass = 2))
    sm$y[sm$y < 0] <- 0
    sm <- sm[-1,]
    avg.sn <- rbind(obs, sm)
    
    
    # Set colors
    col_obs <- "black" 
    col_pred <- "#00BFC4" 
    col_bound <- rgb(col2rgb(col_pred)[1]/255, col2rgb(col_pred)[2]/255, col2rgb(col_pred)[3]/255, 0.4)
    
    # Plot observations for a season
    par(mar = c(5.1, 6.1, 4.1, 2.1), mgp = c(2, 0.5, 0))
    plot(obs.sx, type = "l", lwd = 2, col = col_obs, xlim = c(0, 52), ylim = c(0, maxInc),
         xlab = "Season week", ylab = "", las = 1, tck = -0.04, 
         cex.axis = 1.2, cex.lab = 1.3)
    title(ylab = "Incidence (cases/week)", mgp = c(3, 0.5, 0), cex.lab = 1.3)
    lines(avg.sn$x[(wk+1):52], avg.sn$y[(wk+1):52], type = "l", lwd = 4, col = col_pred)
    polygon(c(upper$x, rev(lower$x) ), c(upper$y, rev(lower$y) ), col = col_bound, border = NA)
    rect(wk, 0, wk+1, maxInc, col = "white", border = NA)
    rect(52, 0, 53, maxInc, col = "white", border = NA)
    lines(obs.sx, type = "l", lwd = 2, col = col_obs)
    abline(v = wk, lwd = 2, lty = 2) # The forecasting week
    
    file_name = paste("../output/", fig,"/forecast_wk_", wk, "_season_", seasonNumber, ".pdf", sep="")
    dev.copy(pdf, file = file_name, width = 4, height = 3)
    dev.off()
  }
}







# Create boxplots for all 7 forecasts for a target and season
boxplot.target <- function(forecast.seasons, seasonNumber, obs.sx, fig, target) {
  
  # Set colors
  col_pred <- "#00BFC4" 
  col_bound <- rgb(col2rgb(col_pred)[1]/255, col2rgb(col_pred)[2]/255, col2rgb(col_pred)[3]/255, 0.4)
  
  sfMax <- 45 # length(previous.simplexes)
  
  ds <- list()
  
  for (i in 1:7) { # For each forecasting week 0, 4, ..., 24
    # Obtain the skew normal distribution based on forecasts
    preds <- forecast.seasons[[seasonNumber]][[i]][1:sfMax]
    
    ds[[i]] <- rsnorm(10000, mean = snormFit(preds)$par[[1]],
                      sd = snormFit(preds)$par[[2]],
                      xi = snormFit(preds)$par[[3]])
    if (target == 1) {
      ds[[i]][ds[[i]] < 1] <- 1
      ds[[i]][ds[[i]] > 52] <- 52
    } else {
      ds[[i]][ds[[i]] < 0] <- 0
    }
  }
  
  # Target-specific parameters
  if (target == 1) {
    targetName <- "peak.wk"
    obs <- which.max(obs.sx)
    y.lim <- c(1, 52)
    y.lab <- "Peak week"
    y.place <- 2
  } else if (target == 2) {
    targetName <- "peak.inc"
    obs <- max(obs.sx)
    y.lim <- c(0, 400)
    y.lab <- "Peak incidence"
    y.place <- 2.8
  } else {
    targetName <- "tot.inc"
    obs <- sum(obs.sx)
    y.lim <- c(0, 8000)
    y.lab <- "Seasonal incidence"
    y.place <- 3.6
  }
  
  par(mar = c(5.1, 6.1, 4.1, 2.1), mgp = c(2, 0.5, 0))
  boxplot(ds, col = col_bound, names = c(0, 4, 8, 12, 16, 20, 24), outline = FALSE, ylim = y.lim,
          xlab = "", ylab = "", las = 1, medcol = col_pred, xaxt = "n", yaxt = "n")
  title(xlab = "Forecasting week", mgp = c(2, 0.5, 0), cex.lab = 1.3)
  title(ylab = y.lab, mgp = c(y.place, 0.5, 0), cex.lab = 1.3)
  axis(1, tck = -0.04, cex.axis = 1.2, cex.lab = 1.3, at = 1:7, labels = c(0, 4, 8, 12, 16, 20, 24))
  axis(2, tck = -0.04, cex.axis = 1.2, cex.lab = 1.3, las = 1)
  abline(h = obs, lwd = 2, lty = 2)
  
  file_name = paste("../output/", fig, "/", targetName, "_season_", seasonNumber, ".pdf", sep="")
  dev.copy(pdf, file = file_name, width = 3.5, height = 3.5)
  dev.off()
}







# END OF SCRIPT