#############################
# Data Formatting Functions #
#                           #
# nicole.nova@stanford.edu  #
#############################


# Load data and remove outlier
data.clean <- function(dataset){
  
  # We load the data and create our main dataframe PR, which stands for "Puerto Rico"
  # since we are studying environmental variables and disease dynamics in San Juan, Puerto Rico
  file_name = paste("../data/", dataset, ".csv", sep="")
  PR <- read.csv(file_name, header = T)
  
  # There is an extreme outlier for precipitation (574.3 mm) that is likely instrumental data error 
  # (over 2x the next highest value) so we will replace it with that next highest value
  prec_max_value <- which.max(PR$prec) # Find the extreme outlier value
  PR$prec[prec_max_value] <- 0 # Remove it
  PR$prec[prec_max_value] <- max(PR$prec) # Replace it with the next highest value
  
  # Then we add a new variable that measures the average weekly temperature in Celsius (°C), 
  # and remove the variable that measures temperature in °F
  PR$temp <- PR$tavg - 273.15
  PR$tavg <- NULL
  
  return(PR)
}





# Add lagged and averaged lagged variables 
data.process <- function(PR){
  
  # Next, we will add a date variable that is formatted for whole time series analyses
  PR$date <- as.Date(as.character(levels(PR$week_start_date)), 
                     format = "%m/%d/%Y")[PR$week_start_date]

  # Clean up order of columns
  # Create an array containing all the variable names
  vbs <- c("cases", "temp", "prec", "mu")
  PR <- PR[, c("date", "week", vbs)]

  # Next, we add new lagged variables
  for (i in 0:25){
    PR[[paste("cases.lag_", i, sep="")]] <- data.table::shift(PR$cases, i)
  }
  
  for (i in 0:25){
    PR[[paste("temp.lag_", i, sep="")]] <- data.table::shift(PR$temp, i)
  }
  
  for (i in 0:25){
    PR[[paste("prec.lag_", i, sep="")]] <- data.table::shift(PR$prec, i)
  }
  
  for (i in 0:25){
    PR[[paste("mu.lag_", i, sep="")]] <- data.table::shift(PR$mu, i)
  }
  
  # Then, we add new averaged lagged variables
  for (i in 0:19){
    PR[[paste("cases.avglag_", i, sep="")]] <- (PR[[paste("cases.lag_", i, sep="")]] +
                                                  PR[[paste("cases.lag_", i+1, sep="")]] +
                                                  PR[[paste("cases.lag_", i+2, sep="")]] +
                                                  PR[[paste("cases.lag_", i+3, sep="")]] +
                                                  PR[[paste("cases.lag_", i+4, sep="")]] +
                                                  PR[[paste("cases.lag_", i+5, sep="")]] +
                                                  PR[[paste("cases.lag_", i+6, sep="")]])/7
  }
  
  for (i in 0:19){
    PR[[paste("temp.avglag_", i, sep="")]] <- (PR[[paste("temp.lag_", i, sep="")]] +
                                                 PR[[paste("temp.lag_", i+1, sep="")]] +
                                                 PR[[paste("temp.lag_", i+2, sep="")]] +
                                                 PR[[paste("temp.lag_", i+3, sep="")]] +
                                                 PR[[paste("temp.lag_", i+4, sep="")]] +
                                                 PR[[paste("temp.lag_", i+5, sep="")]] +
                                                 PR[[paste("temp.lag_", i+6, sep="")]])/7
  }
  
  for (i in 0:19){
    PR[[paste("prec.avglag_", i, sep="")]] <- (PR[[paste("prec.lag_", i, sep="")]] +
                                                 PR[[paste("prec.lag_", i+1, sep="")]] +
                                                 PR[[paste("prec.lag_", i+2, sep="")]] +
                                                 PR[[paste("prec.lag_", i+3, sep="")]] +
                                                 PR[[paste("prec.lag_", i+4, sep="")]] +
                                                 PR[[paste("prec.lag_", i+5, sep="")]] +
                                                 PR[[paste("prec.lag_", i+6, sep="")]])/7
  }
  
  for (i in 0:19){
    PR[[paste("mu.avglag_", i, sep="")]] <- (PR[[paste("mu.lag_", i, sep="")]] +
                                               PR[[paste("mu.lag_", i+1, sep="")]] +
                                               PR[[paste("mu.lag_", i+2, sep="")]] +
                                               PR[[paste("mu.lag_", i+3, sep="")]] +
                                               PR[[paste("mu.lag_", i+4, sep="")]] +
                                               PR[[paste("mu.lag_", i+5, sep="")]] +
                                               PR[[paste("mu.lag_", i+6, sep="")]])/7
  }

  return(PR)
}






# Create a dataframe with seasonal trends
data.season <- function(PR) {
  # The seasonal variables of cases, temp, prec, and mu are averages over all years
  # Create a vector with 52 weeks of a year
  year.week <- c(1:52)
  
  # Create new seasonal (empty) vectors 
  cases.season <- rep(0,52)
  temp.season <- rep(0,52) 
  prec.season <- rep(0,52)
  mu.season <- rep(0,52)
  
  # Calculating year-week averages for all year for cases, temp, prec and mu
  for(i in 1:52){
    data.sub <- subset(PR, week == i)
    cases.season[i] <- mean(data.sub$cases)
    temp.season[i] <- mean(data.sub$temp)
    prec.season[i] <- mean(data.sub$prec)
    mu.season[i] <- mean(na.omit(data.sub$mu))
  }
  
  # Add new variables (log-scaling cases and precipitation due to large value ranges)
  # There are some zeros in the prec data, so we will do log(x + 1)
  cases.season.log <- log10(cases.season) # no zeros
  prec.season.log <- log10(prec.season) # no zeros
  
  # Create new seasonal dataframe
  list.season <- list(year.week, cases.season.log, temp.season, prec.season.log, 
                      cases.season, prec.season, mu.season)
  PR.season <- setNames(do.call(data.frame, list.season), 
                        c("week","cases.season.log", "temp.season",
                          "prec.season.log", "cases.season", "prec.season", "mu.season"))
  
  # Add new variables: lagged temperature, precipitation and mu
  for (i in 0:16){
    PR.season[[paste("cases.season.log.lag_", i, sep="")]] <- data.table::shift(PR.season$cases.season.log, i)
  }
  
  for (i in 0:16){
    PR.season[[paste("temp.season.lag_", i, sep="")]] <- data.table::shift(PR.season$temp.season, i)
  }
  
  for (i in 0:16){
    PR.season[[paste("prec.season.log.lag_", i, sep="")]] <- data.table::shift(PR.season$prec.season.log, i)
  }
  
  for (i in 0:17){
    PR.season[[paste("mu.season.lag_", i, sep="")]] <- data.table::shift(PR.season$mu.season, i)
  }
  
  # Add new variables: lagged averaged variables of temp, prec and mu
  for (i in 0:10){
    PR.season[[paste("temp.season.avglag_", i, sep="")]] <- (PR.season[[paste("temp.season.lag_", i, sep="")]] +
                                                               PR.season[[paste("temp.season.lag_", i+1, sep="")]] +
                                                               PR.season[[paste("temp.season.lag_", i+2, sep="")]] +
                                                               PR.season[[paste("temp.season.lag_", i+3, sep="")]] +
                                                               PR.season[[paste("temp.season.lag_", i+4, sep="")]] +
                                                               PR.season[[paste("temp.season.lag_", i+5, sep="")]] +
                                                               PR.season[[paste("temp.season.lag_", i+6, sep="")]])/7
  }
  
  for (i in 0:10){
    PR.season[[paste("prec.season.log.avglag_", i, sep="")]] <- (PR.season[[paste("prec.season.log.lag_", i, sep="")]] +
                                                                   PR.season[[paste("prec.season.log.lag_", i+1, sep="")]] +
                                                                   PR.season[[paste("prec.season.log.lag_", i+2, sep="")]] +
                                                                   PR.season[[paste("prec.season.log.lag_", i+3, sep="")]] +
                                                                   PR.season[[paste("prec.season.log.lag_", i+4, sep="")]] +
                                                                   PR.season[[paste("prec.season.log.lag_", i+5, sep="")]] +
                                                                   PR.season[[paste("prec.season.log.lag_", i+6, sep="")]])/7
  }
  
  for (i in 0:10){
    PR.season[[paste("mu.season.avglag_", i, sep="")]] <- (PR.season[[paste("mu.season.lag_", i, sep="")]] +
                                                             PR.season[[paste("mu.season.lag_", i+1, sep="")]] +
                                                             PR.season[[paste("mu.season.lag_", i+2, sep="")]] +
                                                             PR.season[[paste("mu.season.lag_", i+3, sep="")]] +
                                                             PR.season[[paste("mu.season.lag_", i+4, sep="")]] +
                                                             PR.season[[paste("mu.season.lag_", i+5, sep="")]] +
                                                             PR.season[[paste("mu.season.lag_", i+6, sep="")]])/7
  }
  
  # Get the last 9 data points (rows 44 - 52) from temp.season, prec.season and append in the begginning
  temp.rest <- PR.season$temp.season[44:52]
  temp.other <- PR.season$temp.season.lag_9[10:52]
  temp.season.lag_9.all <- c(temp.rest, temp.other) 
  
  prec.rest <- PR.season$prec.season.log.avglag_0[44:52]
  prec.other <- PR.season$prec.season.log.avglag_3[10:52]
  prec.season.log.avglag_3.all <- c(prec.rest, prec.other)
  
  mu.rest <- PR.season$mu.season[44:52]
  mu.other <- PR.season$mu.season.lag_17[10:52]
  mu.season.lag_17.all <- c(mu.rest, mu.other)
  
  mu.rest <- PR.season$mu.season[44:52]
  mu.other <- PR.season$mu.season.lag_5[10:52]
  mu.season.lag_5.all <- c(mu.rest, mu.other)
  
  # Average the first 9 values
  prec.1 <- prec.season.log.avglag_3.all[1]
  prec.2 <- prec.season.log.avglag_3.all[2]
  prec.3 <- prec.season.log.avglag_3.all[3]
  prec.4 <- prec.season.log.avglag_3.all[4]
  prec.5 <- prec.season.log.avglag_3.all[5]
  prec.6 <- prec.season.log.avglag_3.all[6]
  prec.7 <- prec.season.log.avglag_3.all[7]
  prec.8 <- prec.season.log.avglag_3.all[8]
  prec.9 <- prec.season.log.avglag_3.all[9]
  prec.10 <- prec.season.log.avglag_3.all[10]
  prec.11 <- prec.season.log.avglag_3.all[11]
  prec.12 <- prec.season.log.avglag_3.all[12]
  prec.13 <- prec.season.log.avglag_3.all[13]
  prec.14 <- prec.season.log.avglag_3.all[14]
  prec.15 <- prec.season.log.avglag_3.all[15]
  
  prec.1.avg <- (prec.1 + prec.2 + prec.3 + prec.4 + prec.5 + prec.6 + prec.7)/7
  prec.2.avg <- (prec.2 + prec.3 + prec.4 + prec.5 + prec.6 + prec.7 + prec.8)/7
  prec.3.avg <- (prec.3 + prec.4 + prec.5 + prec.6 + prec.7 + prec.8 + prec.9)/7
  prec.4.avg <- (prec.4 + prec.5 + prec.6 + prec.7 + prec.8 + prec.9 + prec.10)/7
  prec.5.avg <- (prec.5 + prec.6 + prec.7 + + prec.8 + prec.9 + prec.10 + prec.11)/7
  prec.6.avg <- (prec.6 + prec.7 + prec.8 + + prec.9 + prec.10 + prec.11 + prec.12)/7
  prec.7.avg <- (prec.7 + prec.8 + prec.9 + prec.10 + prec.11 + prec.12 + prec.13)/7
  prec.8.avg <- (prec.8 + prec.9 + prec.10 + prec.11 + prec.12 + prec.13 + prec.14)/7
  prec.9.avg <- (prec.9 + prec.10 + prec.11 + prec.12 + prec.13 + prec.14 + prec.15)/7
  
  prec.rest2 <- c(prec.1.avg, prec.2.avg, prec.3.avg, prec.4.avg, prec.5.avg, prec.6.avg, prec.7.avg, prec.8.avg, prec.9.avg)
  prec.other2 <- prec.season.log.avglag_3.all[10:52]
  prec.season.log.avglag_3.all <- c(prec.rest2, prec.other2)
  
  # Add to seasonal dataframe
  PR.season$temp.season.lag_9.all <- temp.season.lag_9.all
  PR.season$prec.season.log.avglag_3.all <- prec.season.log.avglag_3.all
  PR.season$mu.season.lag_5.all <- mu.season.lag_5.all
  
  return(PR.season)
}





# END OF SCRIPT