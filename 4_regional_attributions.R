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
pacman::p_load(config, data.table, plyr)#, arrow)
# library(data.table)
# library("plyr")

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
foldout_name  = paste0(sSEX, "_", sAGE, "_", sCAU, "_iSEN.", iSEN, "/")
foldout       = paste0( FOLD_DATA_OUT, foldout_name);

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
print(ATTR_PERIOD)

temp_groups = strsplit(CONFIG$TEMP_RANGE, ",")[[1]]
n_groups  = length(temp_groups)

n_regions = length(info_region$code)
eea_subregions = unique(info_region$eea_subregion)
countries = unique(info_region$country_code)
dimnames_agg = list(unlist(list(c("Europe"), eea_subregions, countries)))
ndims_agg = c(1, length(eea_subregions), length(countries)) 
nregions_agg = sum(ndims_agg)

n_simu    = CONFIG$N_SIMU ; if( n_simu < 1000 ){ print( paste0( "  WARNING: Only ", n_simu, " simulations used for AN estimation" ) ); }
col_names = c( sprintf("simu_%s", seq(1:n_simu)), list("attr") )

conf_int = c(0.025,0.975); # two-tailed 95% Confidence Interval # TO DO move to config


# ------------------------------------------------------------------------------
print("Attributions for prediction period adjusted by mortality -")
# ------------------------------------------------------------------------------

if( "Whole Period" %in% ATTR_PERIOD ){
  
  periods = c("Whole Period")
  print( paste0("Analysis for - ", periods) )
  n_periods = length(periods);
  attr_num = attr_frac = array( NA, dim   =   c( n_regions,        n_groups,    1 + n_simu),
                                dimnames = list( info_region$code, temp_groups, col_names ) )
  
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

  # summaries and confidence intervals by regions and countries
  # -----------------------------------------------------------  
  attr_num_agg = attr_frac_agg = array( NA, dim  =    c( nregions_agg,      n_groups,    1 + n_simu),
                                        dimnames = list( dimnames_agg[[1]], temp_groups, col_names ) )
  AN_adj_agg = array( 0, dim=nregions_agg, dimnames=dimnames_agg )
  
  print("Attributions for prediction period for each aggregated region")
  for (r_agg in dimnames_agg[[1]]) {
    print(paste0("  ", r_agg))
    if (r_agg == "Europe") {
      attr_num_agg[r_agg,,] = apply(attr_num, c(2,3), sum)
      AN_adj_agg[r_agg] = sum(AN_adj_factor)
      attr_frac_agg[r_agg,,] = (attr_num_agg[r_agg,,] /  AN_adj_agg[r_agg])*100
    } else {
      
      if (r_agg %in% eea_subregions) {
        agg_idx = (info_region$eea_subregion==r_agg)
      } else if (r_agg %in% countries) {
        agg_idx = (info_region$country_code==r_agg)
      } else { print(paste0("Aggregated Region name not in data: ", r_agg))}
      
      if (sum(agg_idx) > 1) {
        attr_num_agg[r_agg,,] = apply(attr_num[agg_idx,,], c(2,3), sum)
      } else if (sum(agg_idx) == 1) {
        attr_num_agg[r_agg,,] = attr_num[agg_idx,,]
      } else { print("No region selected for aggregation.")} 
      AN_adj_agg[r_agg] = sum(AN_adj_factor[agg_idx])
      attr_frac_agg[r_agg,,] = (attr_num_agg[r_agg,,] /  AN_adj_agg[r_agg])*100
      
    }
  }

  # combine all and add confidence interval of estimates
  # -----------------------------------------------------------
  attr_num_all = abind::abind(attr_num_agg, attr_num, along=1)
  AN_adj_factor_all = abind::abind(AN_adj_agg, AN_adj_factor, along=1)
  attr_frac_all = abind::abind(attr_frac_agg, attr_frac, along=1)

  saveRDS(attr_num_all, paste0(foldout, "AN_simu_wholeperiod", ".rds"))
  saveRDS(AN_adj_factor_all, paste0(foldout, "AN_adj_factor_wholeperiod", ".rds"))
  saveRDS(attr_frac_all, paste0(foldout, "AF_simu_wholeperiod", ".rds"))
  
  AN = get_geo_groups(attr_num_all, 0, conf_int, n_groups, save=TRUE, filename="AN_confint_wholeperiod")
  AF = get_geo_groups(attr_frac_all, 2, conf_int, n_groups, save=TRUE, filename="AF_confint_wholeperiod")
  print("All attributions combined.")

} else { 
  print("not required") 
}

