#########################################################
# Functions for Convergent Cross-Mapping (CCM) Analyses #
#                                                       #
# nicole.nova@stanford.edu                              #
#########################################################

# Obtain seasonal or ebisuzaki surrogate time series for null distributions
getCCMresults <- function(nullMethod, Tperiod){
  # We first create surrogates, run CCM with the surrogates, and collect null distributions
  
  # Create surrogates
  
  # For Figure S7
  # driver xmap cases surrogates
  cases.sur <- make_surrogate_data(PR$cases, method = nullMethod, num_surr = num_sur, T_period = Tperiod)
  
  # Add the cases surrogates to PR
  for (j in 1:num_sur){
    PR[[paste("cases.sur_", j, sep="")]] <- cases.sur[,j]
  }
  
  # For Figure 3 or S8
  # cases xmap driver surrogates
  for (i in 2:4) {
    driver.sur <- make_surrogate_data(PR[[i]], method = nullMethod, num_surr = num_sur, T_period = 52)
    assign(paste(names(PR[i]), ".sur", sep=""), driver.sur)
    
    # Add the driver surrogates to PR
    for (j in 1:num_sur){
      PR[[paste(names(PR[i]), ".sur_", j, sep="")]] <- driver.sur[,j]
    }
  }
  
  
  # CCM with null distributions
  
  # For Figure S7
  # Run CCM with each driver and surrogate cases
  # Run CCM for multiple surrogates and bootstraps: driver xmap cases
  for (k in 2:4) {
    # Perform CCM for each driver xmap cases
    for(j in 1:num_sur){
      # Perform CCM for each surrogate
      c <- ccm(PR, E = get(paste("E_", names(PR[k]), sep="")), lib_column = names(PR[k]), 
               target_column = paste("cases.sur_", j, sep=""), lib_sizes = libs, RNGseed = 2301,
               num_samples = bootrun,
               random_libs = T, replace = F,
               tau = Tau, exclusion_radius = excl_rad, tp = TP)
      
      assign(paste("ccmb.", names(PR[k]), "_xmap_cases.sur_", j, sep=""), c, envir = .GlobalEnv)
    }
  }
  
  
  # For Figure 3 or S8
  # Run CCM with cases and each surrogate driver
  # Run CCM for multiple surrogates and bootstraps: cases xmap driver
  for (k in 2:4) {
    # Perform CCM for each cases xmap driver
    for(j in 1:num_sur){
      # Perform CCM for each surrogate
      c <- ccm(PR, E = E_cases, lib_column = "cases", 
               target_column = paste(names(PR[k]), ".sur_", j, sep=""), lib_sizes = libs, 
               RNGseed = 2301, 
               num_samples = bootrun,
               random_libs = T, replace = F, 
               tau = Tau, exclusion_radius = excl_rad, tp = TP)
      
      assign(paste("ccmb.cases_xmap_", names(PR[k]), ".sur_", j, sep=""), c, envir = .GlobalEnv)
    }
  }
  
  save.image(file=paste("../cache/CCM_null_", nullMethod, ".RData", sep=""))
  
  
  # Compute percentiles for seasonal null distributions
  
  # For Figure S7
  for (k in 2:4) {
    # Collect CCM results of all surrogates to obtain null distribution: cases
    null_cases = NULL
    for (j in 1:num_sur){
      local <- eval(parse(text = paste("ccmb.", names(PR[k]), "_xmap_cases.sur_", j, sep="")))
      null_cases <- rbind(null_cases, local)
    }
    
    # Sort by time series library size
    null_cases <- null_cases[order(null_cases$lib_size),]
    
    # Obtain 2.5%, 50% and 97.5% percentiles from null distribution  
    null_cases_q025 <- NULL
    null_cases_q500 <- NULL
    null_cases_q975 <- NULL
    for (i in 1:length(libs)){
      n <- subset(null_cases, lib_size == libs[i])    # quantile(null_cases[,i], 0.025)
      q025 <- quantile(n$rho, 0.025)
      q500 <- quantile(n$rho, 0.500)
      q975 <- quantile(n$rho, 0.975)
      null_cases_q025 <- rbind(null_cases_q025, q025)
      null_cases_q500 <- rbind(null_cases_q500, q500)
      null_cases_q975 <- rbind(null_cases_q975, q975)
    }
    assign(paste("null.", names(PR[k]), "_xmap_cases.q025", sep=""), null_cases_q025, envir = .GlobalEnv)
    assign(paste("null.", names(PR[k]), "_xmap_cases.q500", sep=""), null_cases_q500, envir = .GlobalEnv)
    assign(paste("null.", names(PR[k]), "_xmap_cases.q975", sep=""), null_cases_q975, envir = .GlobalEnv)
  }
  
  
  # For Figure 3 or S8
  for (k in 2:4) {
    # Collect CCM results of all surrogates to obtain null distribution: driver
    null_driver = NULL
    for (j in 1:num_sur){
      local <- eval(parse(text = paste("ccmb.cases_xmap_", names(PR[k]), ".sur_", j, sep="")))
      null_driver <- rbind(null_driver, local)
    }
    
    # Sort by time series library size
    null_driver <- null_driver[order(null_driver$lib_size),]
    
    # Obtain 2.5%, 50% and 97.5% percentiles from null distribution  
    null_driver_q025 <- NULL
    null_driver_q500 <- NULL
    null_driver_q975 <- NULL
    for (i in 1:length(libs)){
      n <- subset(null_driver, lib_size == libs[i])  # quantile(null_driver[,i], 0.025)
      q025 <- quantile(n$rho, 0.025)
      q500 <- quantile(n$rho, 0.500)
      q975 <- quantile(n$rho, 0.975)
      null_driver_q025 <- rbind(null_driver_q025, q025)
      null_driver_q500 <- rbind(null_driver_q500, q500)
      null_driver_q975 <- rbind(null_driver_q975, q975)
    }
    assign(paste("null.cases_xmap_", names(PR[k]), ".q025", sep=""), null_driver_q025, envir = .GlobalEnv)
    assign(paste("null.cases_xmap_", names(PR[k]), ".q500", sep=""), null_driver_q500, envir = .GlobalEnv)
    assign(paste("null.cases_xmap_", names(PR[k]), ".q975", sep=""), null_driver_q975, envir = .GlobalEnv)
  }
  
  save.image(file=paste("../cache/CCM_null_", nullMethod, ".RData", sep=""))
}




# Plot the subplots for CCM figures 3 and S8
plot.ccm <- function(fig){
  
  var.col <- c("black", "red", "blue", "purple")
  
  for (i in 2:4) {
    cases_driver_q <- get(paste("cases_", names(PR[i]), "_q", sep=""))
    null.cases_xmap_driver.q025 <- get(paste("null.cases_xmap_", names(PR[i]), ".q025", sep=""))
    null.cases_xmap_driver.q500 <- get(paste("null.cases_xmap_", names(PR[i]), ".q500", sep=""))
    null.cases_xmap_driver.q975 <- get(paste("null.cases_xmap_", names(PR[i]), ".q975", sep=""))
    
    # Plot forecast skill vs library size 
    # Plot driver xmap cases 
    plot(cases_driver_q[,2] ~ libs, col = var.col[i], type="l", lty = 1, ylim = c(0,1), lwd = 5, 
         cex.axis = 1.8, cex.lab = 1.1, cex.sub = 1.3,
         xlab="", 
         ylab="")
    #xlab="Number of data points included in cross-mapping", 
    #ylab=expression(paste("Correlation coefficient (",rho,")"))) 
    # median predictive skill vs library size 
    
    
    # Null model
    # Plot cases xmap driver surrogates (95% CI)
    polygon(c(libs, rev(libs)), c(null.cases_xmap_driver.q975, rev(null.cases_xmap_driver.q025)), 
            col = rgb(115/255, 115/255, 115/255, 0.4), border = NA)
    lines(null.cases_xmap_driver.q500 ~ libs, col = rgb(115/255, 115/255, 115/255, 0.4), 
          lwd = 5, lty = 1) # median
    
    # Plot cases cross-mapping driver
    polygon(c(libs, rev(libs)), c(cases_driver_q[,1], rev(cases_driver_q[,3])), 
            col = adjustcolor(var.col[i], alpha.f = 0.4), border = NA)
    lines(cases_driver_q[,2] ~ libs, col = var.col[i], lwd = 5, lty = 1) # median on top
    
    print(get(paste(names(PR[i]), "_cases_MK", sep=""))[[2]]) # Nonsensical (t<0)
    print(get(paste("cases_", names(PR[i]), "_MK", sep=""))[[2]]) # Sensical (t>0)
    
    nam <- c("incidence", "temperature", "rainfall", "susceptibles")
    
    legend(50, 1.06, cex = 1.7, pt.cex = 1, bty = "n",
           c(paste(nam[i], " drives incidence?", sep=""),
             paste(nam[i], " null model", sep="")), 
           lty=c(1,1), lwd = c(5,5), col=c(var.col[i], rgb(115/255, 115/255, 115/255, 0.8)))
    
    file_name = paste("../output/", fig,"/ccm_cases_", names(PR[i]), ".pdf", sep="")
    dev.copy(pdf, file=file_name)
    dev.off()
  }
}




# Plot the subplots for CCM Figure S7
plot.ccm.nonsense <- function(fig){
  for (i in 2:4) {
    driver_cases_q <- get(paste(names(PR[i]), "_cases_q", sep=""))
    null.driver_xmap_cases.q025 <- get(paste("null.", names(PR[i]), "_xmap_cases.q025", sep=""))
    null.driver_xmap_cases.q500 <- get(paste("null.", names(PR[i]), "_xmap_cases.q500", sep=""))
    null.driver_xmap_cases.q975 <- get(paste("null.", names(PR[i]), "_xmap_cases.q975", sep=""))
    
    # Plot forecast skill vs library size 
    # Plot driver cross-mapping cases 
    plot(driver_cases_q[,2] ~ libs, col = "black", type = "l", lty = 1, ylim = c(0,1), lwd = 5, 
         cex.axis = 1.8, cex.lab = 1.1, cex.sub = 1.3,
         xlab="", 
         ylab="")
    #xlab="Number of data points included in cross-mapping", 
    #ylab=expression(paste("Cross-map skill (",rho,")"))) 
    # median predictive skill vs library size
    
    # Null model
    polygon(c(libs, rev(libs)), c(null.driver_xmap_cases.q975, rev(null.driver_xmap_cases.q025)), 
            col = rgb(115/255, 115/255, 115/255, 0.4), border = NA)
    lines(null.driver_xmap_cases.q500 ~ libs, col = rgb(115/255, 115/255, 115/255, 0.4), 
          lwd = 5, lty = 1) # median
    
    # Plot cases cross-mapping driver
    polygon(c(libs, rev(libs)), c(driver_cases_q[,1], rev(driver_cases_q[,3])), 
            col = rgb(0, 0, 0, 0.4), border = NA)
    lines(driver_cases_q[,2] ~ libs, col = "black", lwd = 5, lty = 1) # median on top
    
    print(get(paste(names(PR[i]), "_cases_MK", sep=""))[[2]]) # Nonsensical (t<0)
    print(get(paste("cases_", names(PR[i]), "_MK", sep=""))[[2]]) # Sensical (t>0)
    
    nam <- c("incidence", "temperature", "rainfall", "susceptibles")
    
    legend(90, 1.06, cex = 1.7, pt.cex = 1, bty = "n",
           c(paste("incidence drives ", nam[i],"?", sep=""), 
             paste(nam[i], " null model", sep="")), 
           lty=c(1,1), lwd = c(5,5), col=c("black", rgb(115/255, 115/255, 115/255, 0.8)))
    
    file_name = paste("../output/", fig,"/ccm_cases_", names(PR[i]), ".pdf", sep="") 
    dev.copy(pdf, file=file_name)
    dev.off()
  }
}





# Perform Kolmogorov-Smirnov tests between CCM of observations and CCM of null
ks.results <- function(driverName) {
  # Collect CCM results of all surrogates to obtain null distribution
  null_driver = NULL
  for (j in 1:num_sur){
    local <- eval(parse(text = paste("ccmb.cases_xmap_", driverName, ".sur_", j, sep="")))
    null_driver <- rbind(null_driver, local)
  }
  
  # Obtain distribution at 100 time series length for null
  n <- subset(null_driver, lib_size == 100)
  null.ds <- n$rho

  # Get the distribution at 100 time series length for cases_xmap_driver
  cases_xmap_driver <- eval(parse(text = paste("cases_xmap_", driverName, sep="")))
  d <- subset(cases_xmap_driver, lib_size == 100)
  driver.ds <- d$rho
  
  # KS tests
  r <- stats::ks.test(driver.ds, null.ds) # Compare distributions
  p_val <- round(r$p.value, 10)  # P values of KS tests
  print(p_val)
  
  d <- round(r$statistic, 3)   # D statistics of KS tests
  print(d)
}






# END OF SCRIPT