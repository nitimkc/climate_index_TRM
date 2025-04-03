# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra (NM)
# - Two-stage DLNM with mixmeta
# - adapted from 
# -   Joan Ballester as of October 2024
# -   Tomás Janos as of February 25 2025
# -   NM as of March 3 2025 
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
pacman::p_load(config, dlnm, mixmeta, stargazer)

source("exposure_lag_response.R")         # dlnm wrappers
source("FWALD.R");                        # Wald Test of the meta-predictors

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
lag_knots  = logknots( c( MIN_LAG, MAX_LAG ), 3 ) 

DF_SEAS = CONFIG$DF_SEASONALITY # df/yr for the seasonal and long-term trends
if( DF_SEAS <= 0 ) { stop(paste0("Invalid degree of freedom for seasonality ", DF_SEAS)); }


# ------------------------------------------------------------------------------
# OBTAIN LOCATION-SPECIFIC ASSOCIATIONS
# ------------------------------------------------------------------------------

# //NM four?? coefficients for each regions to be used for meta analysis
# //NM covariances to incorporate the uncertainties of the model

n_regions   = length(info_region$code)
region_list = list(info_region$code)
qaic = MMT  =  array( NA, dim=n_regions, dimnames=region_list )
crpred_premeta = vector( "list", n_regions ); names(crpred_premeta) = info_region$code

n_curves = 1 + length(lrf_centile_idx) # Cum. ERF and LRF at min/max
n_curves_names = c("CUMU", sprintf("P%03g", 100*centiles[lrf_centile_idx]))
coeff = covar = vector( "list", n_curves ); names(coeff) = names(covar) = n_curves_names
for (i in 1:n_curves) {
  covar[[i]] = vector( "list", n_regions); names(covar[[i]]) = info_region$code
}

# for each curve in each region number of coeffs
if (SPLINE_TYPE == "ns") { n_coeff_ERF = length(temp_knots) + 1
} else                   { n_coeff_ERF = length(temp_knots) + spline_degree }
n_coeff_LRF = length(lag_knots) + 2 
for (i in 1:n_curves) {
  if (i == 1) { coeff[[i]] = matrix(NA, n_regions, n_coeff_ERF, dimnames=region_list)
  } else      { coeff[[i]] = matrix(NA, n_regions, n_coeff_LRF, dimnames=region_list)
  }
}

print("Calculating Location-Specific Associations")
# -------------------------------------------------
temp_quantiles =  lapply(mort_temp_calib, function(x) quantile( x$temp, centiles, na.rm=TRUE ))

start_time = Sys.time() # 48.14179 mins remote
for (r in 1:n_regions) {
  reg = info_region$code[r]
  print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", reg, ")" ) );
  data_calib = calib_mort_temp[[r]]
  
  # 1. models, predictions, cross-basis, error
  model_seasonality = seasonality_model(reg, data_calib, DF_SEAS)
  model_cbasis      = crossbasis_model(data_calib, DF_SEAS, temp_knots, SPLINE_TYPE, 
                                       spline_degree, MIN_LAG, MAX_LAG, lag_knots)
  calib_mort_temp[[r]]$pred_seas   = model_seasonality$prediction
  calib_mort_temp[[r]]$pred_cbasis = model_cbasis$prediction
  
  cb_temp = model_cbasis$crossbasis
  qaic[r] = model_cbasis$QAIC 
  
  # 2. premeta cross predictions, MMT 
  temp_quantiles  = quantile( data_calib$temp, centiles, na.rm=TRUE )
  cross_pred = crosspred_premeta(reg, cb_temp, model_cbasis$model, 
                                 all_temp[[r]], temp_quantiles, data_calib$temp, 
                                 CONFIG$LOCAL_MIN_MMT, min_PMMT, max_PMMT)
  MMT[[r]]            = cross_pred$MMT
  crpred_premeta[[r]] = cross_pred$cross_pred
  
  # 3. Reduced Coefficients and Covariance
  creduced_erf     = cross_reduced(cb_temp, model_cbasis$model, MMT[[r]])  # for Cum. ERF
  i = 1
  coeff[[1]][r, ]  = creduced_erf$coeff # ?? To DO if $ used instead of numeric index does it work?
  covar[[1]][[r]]  = creduced_erf$vcov
  
  for( i in 2:n_curves ) {                                     # for each LRF # length(lrf_centile_idx)
    predictor_quantile  = quantile( data_calib$temp, centiles[lrf_centile_idx[i]], na.rm = TRUE ) 
    predictor_value     = predictor_quantile + ( (-1) ^ ( centiles[lrf_centile_idx[i]] < 0.500 ) ) * 0.1 * (predictor_quantile == MMT[[r]])
    creduced_lrf        = cross_reduced(cb_temp, model_cbasis$model, MMT[[r]], type="var", value=predictor_value)
    coeff[[i]][r,]  = creduced_lrf$coeff 
    covar[[i]][[r]] = creduced_lrf$vcov 
  }
}
end_time = Sys.time()
print(end_time-start_time)

print( paste0( "Akaike's Information Criterion for Overdispersed Count Data: ", round(sum(qaic)) ) )
cb_temp_arglag = attr( model_cbasis$crossbasis, "arglag" )

# save results
saveRDS( calib_mort_temp, paste0(foldout, "calib_mort_temp.rds") ) # save predictions
saveRDS( qaic,            paste0(foldout, "qaic.rds") )
saveRDS( MMT,             paste0(foldout, "MMT.rds") )
saveRDS( crpred_premeta,  paste0(foldout, "crpred_premeta.rds") )
saveRDS( coeff,           paste0(foldout, "coeff.rds") )
saveRDS( covar,           paste0(foldout, "covar.rds") )
saveRDS( cb_temp_arglag,  paste0(foldout, "cb_temp_arglag.rds") )
# saveRDS( creduced_erf,   paste0(foldout, "creduced_erf.rds") )
# saveRDS( creduced_lrf,   paste0(foldout, "creduced_lrf.rds") )


# ------------------------------------------------------------------------------
# BLUP for Cum. ERF and LRF at min/max
# ------------------------------------------------------------------------------

# meta predictors
temp_ann = sapply( mort_temp_calib, function(x) mean(x$temp, na.rm=TRUE) ) # annual temp 
temp_IQR = sapply( mort_temp_calib, function(x)  IQR(x$temp, na.rm=TRUE) ) # temp quantile

# temp_win = sapply( mort_temp_calib, function(x) mean( x$temp[12<=month(x$date) | month(x$date)<=2], na.rm=TRUE ) ); # ??
# temp_sum = sapply( mort_temp_calib, function(x) mean( x$temp[ 6<=month(x$date) & month(x$date)<=8], na.rm=TRUE ) );
# temp_p01 = sapply( mort_temp_calib, function(x) quantile( x$temp, 0.01, na.rm=TRUE ) ); 
# temp_p99 = sapply( mort_temp_calib, function(x) quantile( x$temp, 0.99, na.rm=TRUE ) ); 
# names(temp_p01) = names(temp_p99) = info_region$code
# ?? TO DO - get correlation between metapredictors

# mixmeta and BLUPs for each curve
mvar_postmeta = blup_postmeta = vector("list", n_curves)
names(mvar_postmeta) = names(blup_postmeta) = names(coeff) 

start_time = Sys.time() # 1.326498 hours, 33.07614 mins remote
for( i in 1:length(coeff) ) {
  print( paste0("Multivariate Meta-Analysis - ", i, ".", names(coeff)[i]) )
  
  # 1. arguments for mixmeta model
  meta_arg = list( formula = coeff[[i]] ~ temp_ann + temp_IQR, 
                   S       = covar[[i]],
                   data    = data.table(reg=info_region$code),
                   control = list(showiter=TRUE, igls.inititer=10, maxiter=CONFIG$MAX_ITER),
                   method  = "reml" )
  if( as.logical(CONFIG$COUNTRY_RANDOM_EFFECT) ) { 
    meta_arg = c(meta_arg, 
                 list( random  =~ 1 | factor(info_region$country_code)/factor(info_region$code) )) }
  
  # 2. mixmeta model
  mixmeta_model = do.call( mixmeta, meta_arg )
  mvar_postmeta[[i]] = mixmeta_model
  
  model_summary = summary( mixmeta_model)
  # print(model_summary)
  
  mixmeta_coeff = model_summary$lab$p
  print( paste0("Mixmeta coeffs", model_summary$coefficients ) )
  print( paste0("Mixmeta AIC & BIC ", model_summary$AIC, model_summary$BIC ) )
  
  # 3. wald test of the predictors of mixmeta model
  if( length(mixmeta_coeff)  > 1 ) { # ?? why would number of coefficients be less than 1? in fact why would not be the same as in our formula?
    for( p in 1:length(mixmeta_coeff)  ) {
      fwald_estimate = FWALD(mvar_postmeta[[i]], mixmeta_coeff[p])
      print(paste0("  Wald Test of ", mixmeta_coeff[p], ": p = ", 
                   sprintf( "%.10f", fwald_estimate )) ); }
  } else { print(paste0("Mixmeta has only one predictor", mixmeta_coeff)) }
  
  # 4. BLUP for each curve using the mixmeta mvar?? and covar # ??
  blup_postmeta[[i]] = blup( mvar_postmeta[[i]], vcov = TRUE )
}
end_time = Sys.time();
print(end_time-start_time);

# save results
saveRDS( mvar_postmeta, paste0(foldout, "mvar_postmeta.rds") )
saveRDS( blup_postmeta, paste0(foldout, "blup_postmeta.rds") )
# saveRDS( temp_ann,      paste0(foldout, "temp_ann.rds") )
# saveRDS( temp_IQR,      paste0(foldout, "temp_IQR.rds") )
