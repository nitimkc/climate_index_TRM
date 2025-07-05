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
pacman::p_load(config)
# install.packages('DistributionUtils', repos = c('https://dsco036.r-universe.dev', 'https://cloud.r-project.org'))

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

n_curves = 1 + length(lrf_centile_idx) # Cum. ERF and LRF at min/max
n_curves_names = c("CUMU", sprintf("P%03g", 100*centiles[lrf_centile_idx]))

temp_quantiles = lapply(calib_mort_temp, function(x) quantile( x$temp, centiles, na.rm=TRUE ))


# ------------------------------------------------------------------------------
print("Cumulative Exposure and Lag Responses for each region");
# ------------------------------------------------------------------------------

n_regions = length(info_region$code)
MMT_postmeta = array( NA, dim=n_regions, dimnames=list(info_region$code) )
crpred_postmeta = vector( "list", n_curves ); names(crpred_postmeta) = n_curves_names
for( i in 1:n_curves ){
  crpred_postmeta[[i]] = vector("list", n_regions); names(crpred_postmeta[[i]]) = info_region$code;
}

start_time = Sys.time()
for( r in 1:n_regions ){
  reg = info_region$code[r]
  print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", reg, ")" ) );
  data_calib    = calib_mort_temp[[reg]]
  temp_centiles = temp_quantiles[[reg]]
  all_temp_reg = all_temp[[reg]]
  
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
                             at  = all_temp[[reg]],
                             cen = mean(temp_centiles, na.rm=TRUE) )
  if( length(mort_crosspred$predvar) != length(all_temp[[reg]]) ){
    stop( paste0("Length of data and cross prediction (without centering) is not equal for ", r_code, " region !!!" ) ); 
  }
  
  min_MMT =      which( all_temp_reg >= quantile( data_calib$temp, min_PMMT, na.rm=TRUE) )  [1]
  max_MMT = rev( which( all_temp_reg <= quantile( data_calib$temp, max_PMMT, na.rm=TRUE) ) )[1]
  MMT_idx = min_mort_temp(mort_crosspred$allRRfit, CONFIG$LOCAL_MIN_MMT, min_MMT, max_MMT)
  MMT_postmeta[reg] = mort_crosspred$predvar[MMT_idx]
  # MMT_postmeta[reg] = min_mort_temp( mort_crosspred, CONFIG$LOCAL_MIN_MMT, data_calib$temp, 
  #                                  c(min_PMMT,max_PMMT), all_temp[[reg]] )
  
  # with centering for all curves - Cum ERF + LRF at min/max
  print("  postmeta cross-prediction")
  for( i in 1:n_curves ){
    print( paste0( "     ", i, ".  ", names(crpred_postmeta)[i]) )
    if( i == 1 ){ ob = ob_temp ; at = temp_centiles
    } else      { ob = ob_lag  ; at = seq(MIN_LAG, MAX_LAG, 0.1) }
    
    crpred_postmeta[[i]][[reg]] = crosspred(basis = ob,
                                            coef  = blup_postmeta[[i]][[r]]$blup, 
                                            vcov  = blup_postmeta[[i]][[r]]$vcov,
                                            model.link = "log",
                                            at  = at,
                                            cen = MMT_postmeta[reg] ) ;
    if( length(crpred_postmeta[[i]][[reg]]$predvar) != length(at) ){
      stop( paste0("Length of data and cross prediction (with centering) is not equal for ", r_code, " region !!!" ) ); 
    }
  }
}
end_time = Sys.time()
print(end_time-start_time)

# save results
saveRDS( MMT_postmeta,    paste0(foldout, "MMT_postmeta.rds") )
saveRDS( crpred_postmeta, paste0(foldout, "crpred_postmeta.rds") )


# ------------------------------------------------------------------------------
print("Cumulative Exposure and Lag Responses for each aggregated region");
# ------------------------------------------------------------------------------

eea_subregions = unique(info_region$eea_subregion)
countries = unique(info_region$country_code)
dimnames_agg = list(unlist(list(c("Europe"), eea_subregions, countries)))
ndims_agg = c(1, length(eea_subregions), length(countries)) 
nregions_agg = sum(ndims_agg)

MMT_postmeta_agg = array( NA, dim=nregions_agg, dimnames=dimnames_agg )  
crpred_postmeta_agg = vector( "list", n_curves ); names(crpred_postmeta_agg) = n_curves_names
for( i in 1:n_curves ){
  crpred_postmeta_agg[[i]] = vector("list", nregions_agg); names(crpred_postmeta_agg[[i]]) = dimnames_agg[[1]];
}

start_time = Sys.time()
# for each regionally aggregated data
for (r_agg in dimnames_agg[[1]]) {
  print(paste0("  ", r_agg))
  if (r_agg == "Europe") {
    temp_centiles = colMeans(do.call(rbind,temp_quantiles))
  } else {
    if (r_agg %in% eea_subregions) {
      agg_idx = (info_region$eea_subregion==r_agg)
    } else if (r_agg %in% countries) {
      agg_idx = (info_region$country_code==r_agg)
    } else {
      print(paste0("Aggregated Region name not in data: ", r_agg)) 
    }
    temp_centiles = colMeans(do.call(rbind,temp_quantiles[agg_idx]))
  }

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
  mort_crosspred = ob_temp %*% blup_postmeta_agg[[i]][[r_agg]]$fit
  
  min_MMT = which(centiles==min_PMMT)
  max_MMT = which(centiles==max_PMMT)
  MMT_idx = min_mort_temp( mort_crosspred, CONFIG$LOCAL_MIN_MMT, min_MMT, max_MMT )
  PMMT = pmin( pmax(100 * centiles[MMT_idx], 0), 100 )  # Percentile of the MMT aggregated regions
  if (DistributionUtils::is.wholenumber(PMMT)) { 
    MMT_postmeta_agg[r_agg] = temp_centiles[ paste0( PMMT, ".0%" ) ]
  } else { 
    MMT_postmeta_agg[r_agg] = temp_centiles[ paste0( PMMT,   "%" ) ]; 
  }
  
  # with centering for all curves - Cum ERF + LRF at min/max
  print("  Postmeta cross-prediction")
  for( i in 1:n_curves ){
    print( paste0( "     ", i, ".  ", names(crpred_postmeta)[i]) )
    if( i == 1 ){ ob = ob_temp ; at = temp_centiles
    } else      { ob = ob_lag  ; at = seq(MIN_LAG, MAX_LAG, 0.1) }
    
    crpred_postmeta_agg[[i]][[r_agg]] = crosspred(basis = ob,
                                                  coef  = blup_postmeta_agg[[i]][[r_agg]]$fit,
                                                  vcov  = blup_postmeta_agg[[i]][[r_agg]]$vcov,
                                                  model.link = "log",
                                                  at  = at,
                                                  cen = MMT_postmeta_agg[r_agg] ) ;
    if( length(crpred_postmeta_agg[[i]][[r_agg]]$predvar) != length(at) ){
      stop( paste0("Length of data and cross prediction (with centering) is not equal for ", r_code, " region !!!" ) ); 
    }
  }
}
end_time = Sys.time()
print(end_time-start_time)

# save results
saveRDS( MMT_postmeta_agg,    paste0(foldout, "MMT_postmeta_agg.rds") )
saveRDS( crpred_postmeta_agg, paste0(foldout, "crpred_postmeta_agg.rds") )
