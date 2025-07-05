# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra (NM)
# - Two-stage DLNM with mixmeta
# - adapted from 
# -   Joan Ballester as of October 2024
# -   Tomás Janos as of February 25 2025
# -   NM as of March 3 2025 
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# REQUIRED LIBRARIES, FUNCTIONS AND CONFIG
# ------------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(config, dlnm, mixmeta, stargazer)

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

n_curves = 1 + length(lrf_centile_idx) # Cum. ERF and LRF at min/max
n_curves_names = c("CUMU", sprintf("P%03g", 100*centiles[lrf_centile_idx]))

temp_quantiles =  lapply(calib_mort_temp, function(x) quantile( x$temp, centiles, na.rm=TRUE ))


# ------------------------------------------------------------------------------
# STAGE 1. OBTAIN LOCATION-SPECIFIC ASSOCIATIONS
# ------------------------------------------------------------------------------

n_regions   = length(info_region$code)
region_list = list(info_region$code)
qaic = MMT  =  array( NA, dim=n_regions, dimnames=region_list )
crpred_premeta = vector( "list", n_regions ); names(crpred_premeta) = info_region$code
coeff = covar = vector( "list", n_curves ); names(coeff) = names(covar) = n_curves_names
for (i in 1:n_curves) {
  covar[[i]] = vector( "list", n_regions); names(covar[[i]]) = info_region$code
}
# for each curve in each region number of coeffs
if (SPLINE_TYPE == "ns") { n_coeff_ERF = length(temp_knots) + 1
} else                   { n_coeff_ERF = length(temp_knots) + spline_degree }
n_coeff_LRF = length(lag_knots) + 2 
for (i in 1:n_curves) {
  if (i == 1) { coeff[[i]] = matrix(NA, n_regions, n_coeff_ERF, dimnames=region_list) # ERF
  } else      { coeff[[i]] = matrix(NA, n_regions, n_coeff_LRF, dimnames=region_list) # LRF
  }
}
# coefficients for each regions to be used for meta analysis
# covariances to incorporate the uncertainties of the model

print("Calculating Location-Specific Associations")
# -------------------------------------------------
start_time = Sys.time() # 48.14179 mins remote
for (r in 1:n_regions) {
  reg = info_region$code[r]
  print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", reg, ")" ) );
  data_calib = calib_mort_temp[[reg]]
  
  # 1. models, predictions, cross-basis, error
  model_seasonality = seasonality_model(reg, data_calib, DF_SEAS)
  model_cbasis      = crossbasis_model(data_calib, DF_SEAS, temp_knots, SPLINE_TYPE, 
                                       spline_degree, MIN_LAG, MAX_LAG, lag_knots)
  calib_mort_temp[[reg]]$pred_seas   = model_seasonality$prediction
  calib_mort_temp[[reg]]$pred_cbasis = model_cbasis$prediction
  
  cb_temp = model_cbasis$crossbasis
  qaic[reg] = model_cbasis$QAIC 
  
  # 2. premeta cross predictions, MMT 
  cross_pred = crosspred_premeta(reg, cb_temp, model_cbasis$model,all_temp[[reg]], 
                                 temp_quantiles[[reg]], data_calib$temp,
                                 CONFIG$LOCAL_MIN_MMT, min_PMMT, max_PMMT)
  MMT[[reg]]            = cross_pred$MMT
  crpred_premeta[[reg]] = cross_pred$cross_pred
  
  # 3. Reduced Coefficients and Covariance
  for( i in 1:n_curves ) {                                     
    if( i == 1 ){ # for Cum. ERF
      creduced     = cross_reduced(cb_temp, model_cbasis$model, MMT[[reg]])
    } else      { # for each LRF # length(lrf_centile_idx)
      predictor_quantile  = quantile( data_calib$temp, centiles[lrf_centile_idx[i-1]], na.rm = TRUE ) 
      predictor_value     = predictor_quantile + ( (-1) ^ ( centiles[lrf_centile_idx[i-1]] < 0.500 ) ) * 0.1 * (predictor_quantile == MMT[[reg]])
      creduced            = cross_reduced(cb_temp, model_cbasis$model, MMT[[reg]], type="var", value=predictor_value)
    }
    coeff[[i]][reg, ]  = creduced$coeff # ?? To DO if $ used instead of numeric index does it work?
    covar[[i]][[reg]]  = creduced$vcov
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


# ------------------------------------------------------------------------------
# STAGE 2. Mixmeta and BLUP for Cum. ERF and LRF at min/max
# ------------------------------------------------------------------------------

# mixmeta and BLUPs for each curve
mvar_postmeta = blup_postmeta = vector("list", n_curves)
names(mvar_postmeta) = names(blup_postmeta) = n_curves_names

# meta predictors
temp_ann = sapply( calib_mort_temp, function(x) mean(x$temp, na.rm=TRUE) ) # annual temp 
temp_IQR = sapply( calib_mort_temp, function(x)  IQR(x$temp, na.rm=TRUE) ) # temp quantile

# temp_win = sapply( calib_mort_temp, function(x) mean( x$temp[12<=month(x$date) | month(x$date)<=2], na.rm=TRUE ) ); # ??
# temp_sum = sapply( calib_mort_temp, function(x) mean( x$temp[ 6<=month(x$date) & month(x$date)<=8], na.rm=TRUE ) );
# temp_p01 = sapply( calib_mort_temp, function(x) quantile( x$temp, 0.01, na.rm=TRUE ) ); 
# temp_p99 = sapply( calib_mort_temp, function(x) quantile( x$temp, 0.99, na.rm=TRUE ) ); 
# names(temp_p01) = names(temp_p99) = info_region$code
# ?? TO DO - get correlation between metapredictors

print("Fitting Mixmeta")
# ----------------------
start_time = Sys.time() # 1.326498 hours, 33.07614 mins remote
for (i in 1:n_curves ) {
  print( paste0("Multivariate Meta-Analysis - ", i, ".", n_curves_names[i]) )
  
  # 1. mixmeta model
  if( as.logical(CONFIG$COUNTRY_RANDOM_EFFECT) ) { 
    mixmeta_model = mixmeta( formula = coeff[[i]] ~ temp_ann + temp_IQR, 
                             S       = covar[[i]],
                             data    = data.table(reg=info_region$code),
                             control = list(showiter=TRUE, igls.inititer=10, maxiter=CONFIG$MAX_ITER),
                             method  = "reml",
                             random =~ 1 | factor(info_region$country_code)/factor(info_region$code) )
  } else {
    mixmeta_model = mixmeta( formula = coeff[[i]] ~ temp_ann + temp_IQR, 
                             S       = covar[[i]],
                             data    = data.table(reg=info_region$code),
                             control = list(showiter=TRUE, igls.inititer=10, maxiter=CONFIG$MAX_ITER),
                             method  = "reml")
  }
  
  # # 1. arguments for mixmeta model
  # meta_arg = list( formula = coeff[[i]] ~ temp_ann + temp_IQR, 
  #                  S       = covar[[i]],
  #                  data    = data.table(reg=info_region$code),
  #                  control = list(showiter=TRUE, igls.inititer=10, maxiter=CONFIG$MAX_ITER),
  #                  method  = "reml" )
  # if( as.logical(CONFIG$COUNTRY_RANDOM_EFFECT) ) { 
  #   meta_arg = c(meta_arg, 
  #                list( random  =~ 1 | factor(info_region$country_code)/factor(info_region$code) )) }
  # 
  # # 2. mixmeta model
  # mixmeta_model = do.call( mixmeta, meta_arg )
  
  mvar_postmeta[[i]] = mixmeta_model
  model_summary = summary( mixmeta_model)
  # print(model_summary)
  
  mixmeta_coeff = model_summary$lab$p
  print( paste0("Mixmeta coeffs", model_summary$coefficients ) )
  print( paste0("Mixmeta AIC & BIC ", model_summary$AIC, model_summary$BIC ) )
  
  # 2. wald test of the predictors of mixmeta model
  if( length(mixmeta_coeff)  > 1 ) { 
    for( p in 1:length(mixmeta_coeff)  ) {
      fwald_estimate = FWALD(mvar_postmeta[[i]], mixmeta_coeff[p])
      print(paste0("  Wald Test of ", mixmeta_coeff[p], ": p = ", 
                   sprintf( "%.10f", fwald_estimate )) ); }
  } else { print(paste0("Mixmeta has only one predictor", mixmeta_coeff)) }
  
  # 3. BLUP for each curve using the mixmeta mvar?? and covar # ??
  blup_postmeta[[i]] = blup( mvar_postmeta[[i]], vcov = TRUE )
  names(blup_postmeta[[i]]) = info_region$code
}
end_time = Sys.time();
print(end_time-start_time);

# save results
saveRDS( mvar_postmeta, paste0(foldout, "mvar_postmeta.rds") )
saveRDS( blup_postmeta, paste0(foldout, "blup_postmeta.rds") )
# saveRDS( temp_ann,      paste0(foldout, "temp_ann.rds") )
# saveRDS( temp_IQR,      paste0(foldout, "temp_IQR.rds") )


# ------------------------------------------------------------------------------
# BLUP for Aggregated Regions (Europe + Countries)
# ------------------------------------------------------------------------------

eea_subregions = unique(info_region$eea_subregion)
countries = unique(info_region$country_code)
dimnames_agg = list(unlist(list(c("Europe"), eea_subregions, countries)))
ndims_agg = c(1, length(eea_subregions), length(countries)) 
nregions_agg = sum(ndims_agg)

blup_postmeta_agg = vector( "list", n_curves ); names(blup_postmeta_agg) = n_curves_names
for (i in 1:n_curves) {
  blup_postmeta_agg[[i]] = vector("list", nregions_agg); names(blup_postmeta_agg[[i]]) = dimnames_agg[[1]]
}

print("Obtaining BLUP for Aggregated Regions using mixmeta")
# ----------------------------------------------------------
start_time = Sys.time()
for (i in 1:n_curves ) {
  print( paste0("Multivariate Meta-Analysis - ", i, ".", n_curves_names[i]) )
  
  # meta predictors for each regionally aggregated data
  for (r_agg in dimnames_agg[[1]]) {
    print(paste0("  ", r_agg))
    if (r_agg == "Europe") {
      temp_ann_agg = mean(temp_ann)
      temp_IQR_agg = mean(temp_IQR)
    } else {
      if (r_agg %in% eea_subregions) {
        agg_idx = (info_region$eea_subregion==r_agg)
      } else if (r_agg %in% countries) {
        agg_idx = (info_region$country_code==r_agg)
      } else {
        print(paste0("Aggregated Region name not in data: ", r_agg)) 
      }
      temp_ann_agg = mean( temp_ann[agg_idx] )
      temp_IQR_agg = mean( temp_IQR[agg_idx] )
    }
    new_data = data.table(temp_ann = mean( temp_ann_agg),
                          temp_IQR = mean( temp_IQR_agg ) )
    
    # 1. get predictions on this new aggregated data using mvar postmeta 
    blup_postmeta_agg[[i]][[r_agg]] = predict(mvar_postmeta[[i]], new_data, vcov=TRUE, format="list")
  }
}
end_time = Sys.time();
print(end_time-start_time);

# save results
saveRDS( blup_postmeta_agg, paste0(foldout, "blup_postmeta_agg.rds") )
