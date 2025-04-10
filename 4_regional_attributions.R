# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra
# - Regional BLUP and Attributables
# - adapted from 
# -   Joan Ballester as of October 2024
# -   TJ as of February 25 2025
# -   NM as of March 5 2025 
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# REQUIRED LIBRARIES, FUNCTIONS AND CONFIG
# ------------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(config, arrow)

source("exposure_lag_response.R")         # dlnm wrappers
source("attributable_morality.R")         # attr wrappers

CONFIG <- config::get() 


# ------------------------------------------------------------------------------
# SET PARAMETERS + READ DATA
# ------------------------------------------------------------------------------

# file paths
FOLD_DATA_IN  = paste0( CONFIG$ROOT, "indata/")
FOLD_DATA_OUT = paste0( CONFIG$ROOT, "outdata/data_out_daily/")

# attributes
sSEX = CONFIG$SEX
sAGE = CONFIG$AGE
sCAU = CONFIG$CAUSE
iSEN = CONFIG$SENSITIVITY # Change when ready for sensitivity analyses

# attribute specific folder
foldout_name  = paste0(sSEX, "_", sAGE, "_", sCAU, "_iSEN.", iSEN)
foldout       = paste0( FOLD_DATA_OUT, foldout_name, "/" );

# load required outputs from previous runs
info_region     = readRDS( paste0(foldout, "info_region.rds") )
calib_mort_temp = readRDS( paste0(foldout, "calib_mort_temp.rds") )
blup_postmeta   = readRDS( paste0(foldout, "blup_postmeta.rds") )
MMT_postmeta    = readRDS( paste0(foldout, "MMT_postmeta.rds") )
pred_mort_temp  = readRDS( paste0(foldout, "pred_mort_temp.rds") )   # mort reqd. for AN
pred_temp_MOV   = readRDS( paste0(foldout, "pred_temp_MOV.rds") )    # AF when mort not avail.

# Three curves required - 
#   1. Exposure-Response Function (ERF) - Cumulative 
#   2. Lag-Response Function (LRF)      - at min percentile temperature
#   3. Lag-Response Function (LRF)      - at max percentile temperature

# min and max centiles
min_PMMT      = CONFIG$MIN_PMMT
max_PMMT      = CONFIG$MAX_PMMT 

# parameters of ERF & LRF
SPLINE_TYPE = CONFIG$TEMP_SPLINE_TYPE
if ( !SPLINE_TYPE %in% c("ns","bs") ) { stop(paste0("Invalid spline function ", SPLINE_TPYE)); }

spline_degree = as.integer(CONFIG$TEMP_SPLINE_DEG);
if ( SPLINE_TYPE=="bs" & is.na(spline_degree) ) {
  stop( paste0("Invalid degree ", spline_degree, "for spline type ", SPLINE_TYPE) ); }

MIN_LAG = CONFIG$LAG_MIN; MAX_LAG = CONFIG$LAG_MAX
if( MIN_LAG <     0    ) { stop(paste0("Invalid value for minimum lag: ", MIN_LAG)); }
if( MAX_LAG <= MIN_LAG ) { stop(paste0("Maximum lag value: ", MAX_LAG, " <= to minimum lag: ", MIN_LAG)); }

temp_knots = strtoi(strsplit(CONFIG$TEMP_KNOTS, ",")[[1]]) / 100

# others
ATTR_PERIOD = strsplit(CONFIG$ATTRIBUTION_PERIOD, ",")[[1]]
temp_groups = strsplit(CONFIG$TEMP_RANGE, ",")[[1]]
conf_int = c(0.025,0.975); # two-tailed 95% Confidence Interval # TO DO move to config


# ------------------------------------------------------------------------------
print("Attributions for prediction period adjusted by mortality")
# ------------------------------------------------------------------------------

if( "Whole Period" %in% ATTR_PERIOD ){
  
  periods = c("Whole Period")
  print( paste0("Analysis for - ", periods) )

  n_regions = length(info_region$code)
  n_periods = length(periods);
  n_groups  = length(temp_groups)
  n_simu    = CONFIG$N_SIMU ; if( n_simu < 1000 ){ print( paste0( "  WARNING: Only ", n_simu, " simulations used for AN estimation" ) ); }
  
  col_names = c( sprintf("simu_%s", seq(1:n_simu)), list("attr") )
  attr_num = attr_frac = array( NA, dim   =   c( n_regions,        n_groups,    1 + n_simu),
                                dimnames = list( info_region$code, temp_groups, col_names ) )
  # AN_CI = AF_CI = array( 0, dim   =    c( n_regions,      n_periods,    n_groups, 2                 ),
  #                        dimnames = list( info_region$code, periods, temp_groups, c("low", "upper") ) )
  
  # Adjustment of missing values (see Function "attrdl.R" and Gasparrini and Leone 2014)
  AN_adj_factor = array( 0, dim=n_regions, dimnames=list(info_region$code) )
  
  pred_thresholds = sapply( pred_mort_temp, function(x) quantile(x$temp, conf_int, na.rm = TRUE) )
  
  start = Sys.time() # 15.16536 mins remote
  for( r in 1:n_regions ){
    print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", info_region$code[r], ")" ) )
    
    data_calib = calib_mort_temp[[r]]
    data_pred  = pred_mort_temp[[r]]
  
    # One-Basis of temperature centered at MMT
    ob_temp_centered = MMT_centered_onebasis(data_calib$temp,
                                             data_pred$temp,
                                             MMT_postmeta[[r]],
                                             temp_knots, SPLINE_TYPE, spline_degree)
    # daily ts of the lagged mortality
    lagged_mort     = tsModel::Lag( data_pred$mort, -seq(MIN_LAG,MAX_LAG) ) # t x n_lags
    avg_lagged_mort = rowMeans( lagged_mort, na.rm = FALSE )                # t x 1 (entire period incl. NAs)
  
    i = 1 # for CUM. ERF
    # TO DO ADD COMMENT
    AN_ts_simu = simulated_attributable_values(ob_temp_centered,
                                               blup_postmeta[[i]][[r]]$blup,
                                               blup_postmeta[[i]][[r]]$vcov,
                                               CONFIG$SEED, n_simu,
                                               CONFIG$LOCAL_MIN_MMT,
                                               min_PMMT, max_PMMT,
                                               col_names, avg_lagged_mort)
    
    time_idx = 1:length( data_pred$date )
    N = length(time_idx)
    if( N <= 0 ){ stop( paste0("No time period selected for region - ", r) ); } # when will this condition exist?
    rownames(AN_ts_simu) = data_pred$date                                       # format(time_idx, "%Y-%m-%d")
  
    # average deaths across lags not counting NAs unlike lagged_mort
    AF_adj_factor    = sum(      avg_lagged_mort[time_idx]                           , na.rm=TRUE );
    AN_adj_factor[r] = sum( rowMeans(lagged_mort[time_idx, , drop=FALSE], na.rm=TRUE), na.rm=TRUE ); # (specific period without NAs)

    # Attributions for each temperature range group
    idx_groups = sapply(temp_groups, group_threshold, data_pred$temp, MMT_postmeta[[r]], pred_thresholds[,r], time_idx)
    if (any(length(idx_groups) <= 0)) { stop("length of index for one or more temperature range group is less than 1.")
    } else {
      # obtain attribution fractions and numbers for each group
      AN_temp_group = lapply(idx_groups, function(x) AN_ts_simu[x,])
      attributions = t(sapply(AN_temp_group, function(x) colSums(x, na.rm=TRUE) / AF_adj_factor))
      stopifnot(
        "Invalid attributions !!!"   = is.finite( attributions )
      )
      # Adjusted Factor (in annual deaths) for Missing Values in the Calculation of the Regional Attributable Number
      # (see Function "attrdl.R" and Gasparrini and Leone 2014):
      AN_adj_factor[r] = 365.25 * AN_adj_factor[r] / N
      
      attr_frac[r, , ] = attributions * 100
      attr_num[r, , ]  = attributions * AN_adj_factor[r]
    }
  }
  end = Sys.time()
  print(end-start)

  # save
  saveRDS( attr_num,      paste0(foldout, "attr_num.rds") )
  saveRDS( attr_frac,     paste0(foldout, "attr_frac.rds") )
  saveRDS( AN_adj_factor, paste0(foldout, "AN_adj_factor.rds") )
  
} else { 
  print("Adjusted attributions for whole period not required") 
  }


# ------------------------------------------------------------------------------
print("AF for prediction period for Modes of Variability (MOV) analysis")
# ------------------------------------------------------------------------------

if('NUTS' %in% ATTR_PERIOD){
  
  # periods = seq( as.Date(CONFIG$NUTS_STARTDATE), as.Date(CONFIG$NUTS_ENDDATE)-MAX_LAG, 1 )
  periods = seq( as.Date(CONFIG$NUTS_STARTDATE), as.Date(CONFIG$NUTS_ENDDATE), 1 )
  periods = format(periods, "%Y-%m-%d")

  n_periods = length(periods)
  n_groups  = length(temp_groups)
  n_simu    = CONFIG$N_SIMU ; if( n_simu < 1000 ){ print( paste0( "  WARNING: Only ", n_simu, " simulations used for AN estimation" ) ); }
  col_names = c( sprintf("simu_%s", seq(1:n_simu)), list("attr") )
  
  pred_thresholds_MOV = sapply( pred_temp_MOV, function(x) quantile( x$temp, conf_int, na.rm = TRUE ) )
  
  # attr_frac = array( NA, dim=c(n_regions,n_simu+1,n_groups), dimnames=list(info_region$code,col_names,temp_groups) )
  # attr_frac too large object to save in memory
  
  # foldout_attr = paste0( foldout, "AF_ts_simu" )
  foldout_attr = paste0( foldout, "AF_ts_simu_pqt" ) # TO DO - MOVE FOLDER NAME TO CONFIG
  if( !file_test( "-d", foldout_attr ) ){ dir.create( file.path(foldout_attr)); }
  
  start = Sys.time() # 15.16536 mins remote
  for( r in 1:n_regions ){
    print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", info_region$code[r], ")" ) )
    
    data_calib = calib_mort_temp[[r]]
    data_pred  = pred_temp_MOV[[r]]
    
    # One-Basis of temperature centered at MMT
    ob_temp_centered = MMT_centered_onebasis(data_calib$temp, 
                                             data_pred$temp,
                                             MMT_postmeta[[r]],
                                             temp_knots, SPLINE_TYPE, spline_degree)
    i = 1 # for CUM. ERF
    # ADD COMMENT
    AF_ts_simu = simulated_attributable_values(ob_temp_centered,
                                               blup_postmeta[[i]][[r]]$blup,
                                               blup_postmeta[[i]][[r]]$vcov,
                                               CONFIG$SEED, n_simu,
                                               CONFIG$LOCAL_MIN_MMT,
                                               min_PMMT, max_PMMT,
                                               col_names) # no mortality adjustment for prediction period
    AF_ts_simu = AF_ts_simu * 100
    rownames(AF_ts_simu) = periods
    
    # Attributions for each temperature range group
    idx_groups = sapply(temp_groups, group_threshold, data_pred$temp, MMT_postmeta[[r]], pred_thresholds_MOV[,r], 1:n_periods)
    filter_cols = array(FALSE, dim=c(n_periods,n_groups), dimnames=list(periods, temp_groups))
    for (g in 1:n_groups) { filter_cols[idx_groups[[g]] , g] = TRUE }
    AF_ts_simu = do.call(dplyr::bind_cols, list(AF_ts_simu, filter_cols)) 
    
    pqt_filepath = tempfile(tmpdir=foldout_attr, fileext=".gzip.parquet")
    write_parquet(AF_ts_simu, sink=pqt_filepath, compression="gzip")           # faster than snappy
    # saveRDS( AF_ts_simu, paste0(foldout_attr, info_region$code[r], ".rds") ) # TO DO current saved versions do not have rownames add if re-run
  }
  rm(filter_cols, AF_ts_simu)
  end = Sys.time()
  print(end-start)

} else {
  print("Attributions for Prediction Period not required")
}

