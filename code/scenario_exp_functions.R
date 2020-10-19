###############################################
# Functions for Scenario Exploration Analyses #
#                                             #
# nicole.nova@stanford.edu                    #
###############################################


# Takes an input matirx or data frame and creates a block of lag coordinates for state space reconstruction
make_block <- function(data, cols, delays, 
                       lib = c(1, NROW(data)), 
                       diff_col = rep(FALSE, length(cols)))
{ 
  # INPUTS:
  #   data - matrix, array, or data.frame with time series variables arranged in columns
  #   cols - vector indices or names of the columns of 'data' to use for each column of block
  #   delays - vector with same length as cols specifying the time displacement for that column
  #   diff_col - vector of logical on whether to apply first differencing to this lag coordinate variable
  #
  # OUTPUT:
  #   block - array with length(cols) of columns, and NROW(data) rows
  
  if(!is.numeric(cols)){
    cols <- as.numeric(factor(levels = colnames(data), x = cols))
  }
  
  lib <- matrix(lib, ncol = 2)
  # data <- as.matrix(data)
  
  ncol <- length(cols)
  nrow <- dim(data)[1]
  block <- as.data.frame(array(NA, dim = c(nrow, ncol)))
  names(block) <- 1:ncol
  
  for (i in 1:ncol)
  {
    I <- 1:nrow
    I_delay <- intersect(I, I + delays[i])
    block[I_delay-delays[i], i] <- data[I_delay, cols[i]]
    if (delays[i] < 0){
      # remove data points that fall at start of lib segments
      block[lib[, 1] - (0:(delays[i] + 1)), i] <- NA
      names(block)[i] <- paste(colnames(data)[cols[i]], '_t-', abs(delays[i]), sep="")  
    } else if (delays[i] > 0) {
      # remove data points that fall at end of lib segments
      block[lib[, 2] - (0:(delays[i] - 1)), i] <- NA
      names(block)[i] <- paste(colnames(data)[cols[i]], '_t+', abs(delays[i]), sep="")  
    } else {
      names(block)[i] <- paste(colnames(data)[cols[i]], '_t', sep="")
    }
    
    if (diff_col[i]){
      block[, i] <- c(NA, diff(block[, i]))
    }
  } # i
  
  
  return(block)
}






scenario_exploration <- function(block,
                                 target_column = 1,
                                 explore_column,
                                 columns = 1:NCOL(block),
                                 delta = 0.05,
                                 lib = c(1,NROW(block)),
                                 pred = c(1,NROW(block)),
                                 ...){
  
  pred <- pred + NROW(block)
  
  block.plus <- block 
  block.plus[, explore_column] = block.plus[, explore_column] + delta/2
  
  block.minus <- block
  block.minus[, explore_column] = block.minus[, explore_column] - delta/2
  
  block.plus <- bind_rows(block, block.plus)
  block.minus <- bind_rows(block, block.minus)
  
  results.plus <- block_lnlp(
    block = block.plus,
    target_column = target_column,
    lib = lib,
    pred = pred,
    stats_only = FALSE,
    short_output = TRUE,
    ...)
  
  results.minus <- block_lnlp(
    block=block.minus,
    target_column = target_column,
    lib = lib,
    pred = pred,
    stats_only = FALSE,
    short_output = TRUE,
    ...)
  
  results.plus <- results.plus[[1]]$model_output %>%
    rename(pred.plus = pred)
  
  results.minus <- results.minus[[1]]$model_output %>%
    rename(pred.minus = pred)
  
  results <- full_join(results.plus,
                       results.minus,
                       by = c("time", "obs")) %>%
    #mutate(delta=pred.plus-pred.minus)            # Only necessary when reporting normalized values
    mutate(delta = (pred.plus-pred.minus)/delta )
  
  return(results)
}





# Normalization functions
norm_rMM <- function(v){
  sqrt((v - min(v, na.rm=TRUE)) / ( max(v, na.rm = TRUE) - min(v, na.rm=TRUE) ))
}


norm_MM <- function(v){
  (v - min(v, na.rm=TRUE)) / ( max(v, na.rm = TRUE) - min(v, na.rm=TRUE) ) 
}  


norm_coeff <- function(v, v_mean, v_sd) {
  (v - v_mean) / v_sd
} 





# END OF SCRIPT