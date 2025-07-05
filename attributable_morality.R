# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra
# - Wrapper functions for attributable fraction and number calculation for
#   various temp. periods, regions and thresholds
# - adapted from 
# -   Joan Ballester as of October 2024
# -   Tomás Janos as of February 25 2025
# -   NM as of March 12 2025 
# ------------------------------------------------------------------------------


if (!require("pacman")) install.packages("pacman")
suppressMessages(  pacman::p_load(config, tidyverse) )
# library(magrittr) # needs to be run every time you start R and want to use %>%
# library(dplyr)    # alternatively, this also loads %>%

# source(".R")         # wrappers


group_threshold <- function(temp_grp, temp, MMT, thresholds, time_idx){
  
  if        (temp_grp =="Total")         {threshold = 1:length(time_idx)
  } else if (temp_grp =="Total Cold")    {threshold = which(temp[time_idx] < MMT            )
  } else if (temp_grp =="Total Heat")    {threshold = which(           MMT < temp[time_idx] )
  } else if (temp_grp =="Moderate Cold") {threshold = which(thresholds[[1]]  < temp[time_idx] & temp[time_idx] < MMT          )
  } else if (temp_grp =="Moderate Heat") {threshold = which(           MMT < temp[time_idx] & temp[time_idx] < thresholds[[2]])
  } else if (temp_grp =="Extreme Cold")  {threshold = which(temp[time_idx] < thresholds[[1]]  )
  } else if (temp_grp =="Extreme Heat")  {threshold = which(thresholds[[2]]  < temp[time_idx] )
  } else                            { stop("ERROR: Invalid Temperature Range !!!"); }
  # print(length(threshold))
  return (threshold)
}


simulated_attributable_values <- function(ob_temp, blup_coeff, blup_covar, seed, 
                                          n_simu, local_min_MMT, min_PMMT, max_PMMT,
                                          col_names, avg_lagged_mort=NA){
  # perturbations of BLUP coeffs (dim: beta x n_simu)
  set.seed(CONFIG$SEED)
  simu_coeff = t( MASS::mvrnorm(n_simu, blup_coeff, blup_covar) ) 
  
  # Cross-predicted mortality 
  reference  = 1 - exp(-ob_temp %*% blup_coeff)
  simulated  = 1 - exp(-ob_temp %*% simu_coeff)
  # check negative values
  if( any( reference < 0 ) ){
    print( paste0( "    WARNING: ", 
                   round( 100 * mean( reference < 0 ), digits = 2 ), 
                   "% of Predicted Attributable Fractions Are Negative" ) );
    if (local_min_MMT & min_PMMT <= 0 & max_PMMT >= 1) {
      print(" Negative Predicted Attributable Fractions are theoretically due to lack of MMT precision")
    }
  }
  
  # Time series of attributable number
  if ( (length(avg_lagged_mort)!=1) & any(!is.na(avg_lagged_mort)) ) {
    reference  = reference * avg_lagged_mort     # value
    simulated  = simulated * avg_lagged_mort     # simulated 
  }
  
  # combine and return
  attributions = cbind(simulated, reference)
  colnames(attributions) = col_names
  # return(attributions)
  return(data.table(attributions)) # ?? faster sapply operation
}


get_geo_groups <- function(data, round_by, conf_int, n_groups, save, filename) {
  
  attr = round(data[,,'attr'], round_by)
  attr_confint = round(apply(data[,,1:n_simu], c(1,2), quantile, conf_int), round_by) 
  attr_confint_all = abind::abind(attr, attr_confint, along=1)  
  dimnames(attr_confint_all)[[1]][1] = "attr"

  df = data.frame( matrix(ncol=1 + dim(attr_confint_all)[3], 
                          nrow=    dim(attr_confint_all)[2] ) )
  colnames(df) = c("Region", dimnames(attr_confint_all)[[3]])
  df$Region = dimnames(attr_confint_all)[[2]]
  
  for (g in 1:n_groups) {
    print(g)
    g_vals = apply(attr_confint_all[,,g], 1, function(x) as.vector(t(x)))
    df[,g+1] = paste0(g_vals[,1], " (", g_vals[,2], ", ", g_vals[,3], ")")
  }
  if (save==TRUE){
    write.csv(df, paste0(foldout, filename, ".csv"), row.names = FALSE)
  }
  return(df)
}


get_confidence_interval <- function(attr){
  
  # attr_confint = array( NA, dim=n_dims, dimnames=dim_names )
  actual = attr[ ,n_simu+1]
  actual = array(actual, c(length(actual),1), list(rownames(attr),"attr"))
  
  confint = t(apply(attr, 1, quantile, conf_int, na.rm=TRUE))
  attr_confint = abind::abind(actual, confint, along=2)
  
  return(attr_confint)
}

