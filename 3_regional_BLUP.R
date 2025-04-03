# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra
# - Regional BLUP and Attributables
# - adapted from 
# -   Joan Ballester as of October 2024
# -   TJ as of February 25 2025
# -   NM as of March 5 2025 
# ------------------------------------------------------------------------------


# start_time = Sys.time();
# end_time = Sys.time();
# print(end_time-start_time);

# rm( list = ls() );
# cat("\014");


# ------------------------------------------------------------------------------
# REQUIRED LIBRARIES, FUNCTIONS AND CONFIG
# ------------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(config)

source("exposure_lag_response.R")         # dlnm wrappers

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
all_temp        = readRDS( paste0(foldout, "all_temp.rds") )
cb_temp_arglag  = readRDS( paste0(foldout, "cb_temp_arglag.rds") )
blup_postmeta   = readRDS( paste0(foldout, "blup_postmeta.rds") )

# Three curves required - 
#   1. Exposure-Response Function (ERF) - Cumulative 
#   2. Lag-Response Function (LRF)      - at min percentile temperature
#   3. Lag-Response Function (LRF)      - at max percentile temperature

# prediction centiles (//NM ensure temp is within the same range for each region) 
centiles = sort( unique( c( seq(  0.0,   1.0, 0.1 ),
                            seq(  1.5,   5.0, 0.5 ),
                            seq(  6.0,  94.0, 1.0 ),
                            seq( 95.0,  98.5, 0.5 ),
                            seq( 99.0, 100.0, 0.1 ) ) / 100 ) ) # # for Cum. ERF ??

# min and max centiles
min_PMMT      = CONFIG$MIN_PMMT
max_PMMT      = CONFIG$MAX_PMMT 
lrf_centile_idx = c( which(centiles == min_PMMT), which(centiles == max_PMMT) ) # for LRF

if( min(centiles)!=0 | max(centiles)!=1 | any(0>centiles | centiles>1) ) {
  stop("Invalid centile Vector for the Predictions of the Cumulative Exposure-Response !!!");}
if( !min_PMMT %in% centiles ) { 
  stop("Invalid minimum centile MMT or centile crosspred !!!"); }
if( !max_PMMT %in% centiles ) { 
  stop("Invalid maximum centile MMT or centile crosspred !!!"); }

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


# ------------------------------------------------------------------------------
print("BLUPs for each region");
# ------------------------------------------------------------------------------
n_regions = length(info_region$code)
MMT_postmeta = array( NA, dim=n_regions, dimnames=list(info_region$code) )

n_curves = 1 + length(lrf_centile_idx) # Cum. ERF and LRF at min/max
n_curves_names = c("CUMU", sprintf("P%03g", 100*centiles[lrf_centile_idx]))
crpred_postmeta = vector( "list", n_curves ); names(crpred_postmeta) = n_curves_names
for( i in 1:n_curves ){
  crpred_postmeta[[i]] = vector("list", n_regions); names(crpred_postmeta[[i]]) = info_region$code;
}

temp_quantiles = lapply(calib_mort_temp, function(x) quantile( x$temp, centiles, na.rm=TRUE ))

start_time = Sys.time()
for( r in 1:n_regions ){
  
  print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", info_region$code[r], ")" ) );
  data_calib    = calib_mort_temp[[r]]
  temp_centiles = temp_quantiles[[r]]
  
  # 1. One-Basis 
  #    of Temperature at centiles
  argvar = list( x     = temp_centiles, 
                 knots = temp_centiles[paste0(100*temp_knots, ".0%")], 
                 fun   = SPLINE_TYPE,
                 Bound = range(temp_centiles, na.rm = TRUE) )
  if( SPLINE_TYPE == "bs" ){ argvar =  c(argvar, list( degree = spline_degree )) }
  ob_temp = do.call( onebasis, argvar)
  
  #    of lags at min/max
  arglag = c(list(x=seq(MIN_LAG,MAX_LAG, 1)), cb_temp_arglag) # confirm this is lag range
  ob_lag = do.call( onebasis, arglag)
  
  # 2. postmeta cross predictions, MMT
  i = 1 # without centering for Cum ERF only 
  mort_crosspred = crosspred(basis = ob_temp,
                             coef  = blup_postmeta[[i]][[r]]$blup, 
                             vcov  = blup_postmeta[[i]][[r]]$vcov,
                             model.link = "log",
                             at  = all_temp[[r]],
                             cen = mean(temp_centiles, na.rm=TRUE) )
  if( length(mort_crosspred$predvar) != length(all_temp[[r]]) ){
    stop( paste0("Length of data and cross prediction (without centering) is not equal for ", r_code, " region !!!" ) ); 
  }
  
  MMT_postmeta[r] = min_mort_temp( mort_crosspred, CONFIG$LOCAL_MIN_MMT, data_calib$temp, 
                                   c(min_PMMT,max_PMMT), all_temp[[r]] )
  
  # with centering for all curves - Cum ERF + LRF at min/max
  print("  postmeta cross-prediction")
  for( i in 1:n_curves ){
    print( paste0( "     ", i, ".  ", names(crpred_postmeta)[i]) )
    if( i == 1 ){ ob = ob_temp ; at = temp_centiles
    } else      { ob = ob_lag  ; at = seq(MIN_LAG, MAX_LAG, 0.1) }
    
    crpred_postmeta[[i]][[r]] = crosspred(basis = ob,
                                          coef  = blup_postmeta[[i]][[r]]$blup, 
                                          vcov  = blup_postmeta[[i]][[r]]$vcov,
                                          model.link = "log",
                                          at  = at,
                                          cen = MMT_postmeta[r] ) ;
    if( length(crpred_postmeta[[i]][[r]]$predvar) != length(at) ){
      stop( paste0("Length of data and cross prediction (with centering) is not equal for ", r_code, " region !!!" ) ); 
    }
  }
}
end_time = Sys.time()
print(end_time-start_time)

# save results
saveRDS( MMT_postmeta,    paste0(foldout, "MMT_postmeta.rds") )
saveRDS( crpred_postmeta, paste0(foldout, "crpred_postmeta.rds") )

