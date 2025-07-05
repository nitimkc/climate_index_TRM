# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra
# - Wrapper function for two-stage DLNM process
# - adapted from 
# -   Joan Ballester as of October 2024
# -   Tomás Janos as of February 25 2025
# -   NM as of March 3 2025 
# ------------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
suppressMessages( pacman::p_load(dlnm, splines) )

model_check <- function(model, reg){
  "model: An R model object
   reg: A string representing the region code for which the model is passed
   
   Various checks for model object."
  
  if( !model$converged ){
    print( paste0("    WARNING: The Model Did Not Converge for ", reg, " region.") );
  }
  if( any( !is.finite( summary(model)[[13]] ) ) ){
    print( paste0("    WARNING: Non-Finite Model Coefficients or Standard Errors model for ", reg, " region.") );
    print( summary(model)[[13]] );
  }
  if( any( !is.finite( coef(model) ) ) ){
    print( paste0("    WARNING: Non-Finite Model Coefficients in model for ", reg, " region.") );
    print( coef(model) );
  }
  if( any( !is.finite( vcov(model) ) ) ){
    print( paste0("    WARNING: Non-Finite Model Covariance Matrix in model for ", reg, " region.") );
    print( vcov(model) );
  }
}


seasonality_model <- function(reg, data, degr_f){
  "reg: A string representing the region code for which the model is passed
   data: A table of data for the corresponding region with date and mortality columns
   degr_f: An integer representing degree of freedom
   
   TO DO - ADD OBJECTIVE"
  
  degr_f = round(degr_f*length(data$date)/365.25)
  model  = glm(formula  = mort ~ ns(date, df=degr_f), 
              family    = quasipoisson, 
              data      = data, 
              na.action = "na.exclude")
  model_check(model, reg)
  prediction = predict(model, type="response")
  
  return( list("model"      = model, 
               "prediction" = prediction ) )
}


crossbasis_model <- function(mort_temp, degr_f, temp_knots, func_temp, 
                             func_degree, min_lag, max_lag, lag_knots) {
  "mort_temp: A table of data for the corresponding region with date, temperature and mortality columns
   degr_f: An integer representing degree of freedom
   temp_knots: 
   func_temp: 
   func_degree:
   min_lag: 
   max_lag: 
   lag_knots:
   
   TO DO - ADD OBJECTIVE"
  
  degr_f = round(degr_f*length(mort_temp$date)/365.25)                    
  
  argvar = list( knots = quantile( mort_temp$temp, temp_knots, na.rm=TRUE ),
                 fun   = func_temp,
                 Bound = range( mort_temp$temp, na.rm=TRUE ) )
  if (func_temp == "bs"){ argvar <- c(argvar, degree=degr_temp) }
  
  cb_temp = crossbasis(mort_temp$temp, 
                       c(min_lag, max_lag),
                       arglag = list( knots = lag_knots),
                       argvar = list( knots = quantile( mort_temp$temp, temp_knots, na.rm=TRUE ),
                                      fun   = func_temp,
                                      Bound = range( mort_temp$temp, na.rm=TRUE ) ) )
  
  model = glm(formula   = mort ~ ns(date, df=degr_f) + dow + cb_temp, 
              family    = quasipoisson, 
              data      = mort_temp, 
              na.action = "na.exclude")
  model_check(model)
  prediction = predict(model, type="response")
  model_qaic = QAIC(model)
  
  return( list("model"       = model, 
               "prediction"  = prediction,
               "crossbasis"  = cb_temp,
               "QAIC"        = model_qaic ) )
}

min_mort_temp <- function(allRRfit, local_min_MMT, min_MMT, max_MMT){
  "RRfit: 
   local_min_MMT: 
   pred: 
   min_MMT:
   max_MMT:
   
   TO DO - ADD OBJECTIVE"
  MMT_minima = 1 + which( diff( sign( diff(allRRfit) ) ) == 2 )
  
  if( local_min_MMT & length(MMT_minima) > 0 ){
    # if at least one local minima exists, MMT is one with lowest relative risk
    MMT = MMT_minima[ which.min( allRRfit[MMT_minima] ) ]
  }else{
    # if no local minima exists, MMT is one with the lowest relative risk within a predefined temperature percentile range
    # // NM ?? any value that has the lowest relative risk within the desired temperature range
    MMT = min_MMT - 1 + which.min( allRRfit[min_MMT:max_MMT] )
  }
  
  return(MMT)
}


# min_mort_temp <- function(cross_pred, local_min_MMT, temp, percentile_MMT, all_temp){
#   "cross_pred:
#    local_min_MMT:
#    percentile_MMT:
#    all_temp:
# 
#    TO DO - ADD OBJECTIVE"
# 
#   MMT_minima = 1 + which( diff( sign( diff(cross_pred$allRRfit) ) ) == 2 )
# 
#   if( local_min_MMT & length(MMT_minima) > 0 ){
#     # if at least one local minima exists, MMT is one with lowest relative risk
#     MMT = cross_pred$predvar[ MMT_minima[ which.min( cross_pred$allRRfit[MMT_minima] ) ] ]
#   }else{
#     # if no local minima exists, MMT is one with the lowest relative risk within a predefined temperature percentile range
#     # // NM ?? any value that has the lowest relative risk within the desired temperature range
#     min_MMT =      which( all_temp >= quantile( temp, percentile_MMT[[1]], na.rm=TRUE ) )  [1]
#     max_MMT = rev( which( all_temp <= quantile( temp, percentile_MMT[[2]], na.rm=TRUE ) ) )[1]
#     MMT = cross_pred$predvar[ min_MMT - 1 + which.min( cross_pred$allRRfit[min_MMT:max_MMT] ) ]
#   }
# 
#   return(MMT)
# }

# 
# model =model_cbasis$model
# alltemp=all_temp[[reg]]
# temp=data_calib$temp
# local_min_MMT=CONFIG$LOCAL_MIN_MMT
crosspred_premeta <- function(reg, cb_temp, model, alltemp, temp_quantiles, temp, 
                              local_min_MMT, min_PMMT, max_PMMT){
  "reg: A string representing the region code for which the model is passed
   cb_temp: 
   model: 
   all_temp:
   temp_quantiles:
   temp:
   local_min_MMT:
   min_PMMT:
   max_PMMT:
   
   TO DO - ADD OBJECTIVE"
  
  # without centering
  # -----------------
  cross_pred = crosspred( cb_temp, model,
                          model.link = "log",
                          at         = alltemp, bylag = 1,
                          cen        = mean(temp, na.rm=TRUE) )
  if( length(cross_pred$predvar) != length(alltemp) ){
    stop( paste0("Length of data and cross prediction (without centering) is not equal for ", reg, " region !!!" ) ); 
  }
  
  # with centering
  # --------------
  min_MMT =      which( alltemp >= quantile( temp, min_PMMT, na.rm=TRUE ) )  [1]
  max_MMT = rev( which( alltemp <= quantile( temp, max_PMMT, na.rm=TRUE ) ) )[1]
  MMT_idx = min_mort_temp( cross_pred$allRRfit, local_min_MMT, min_MMT, max_MMT)
  MMT = cross_pred$predvar[MMT_idx]
  
  # MMT = min_mort_temp( cross_pred, local_min_MMT, temp, c(min_PMMT,max_PMMT), all_temp )
  cross_pred = crosspred( cb_temp, model, model.link = "log", bylag = 1,
                          at         = temp_quantiles,
                          cen        = MMT); # MMT to be calculated at temperature centiles, rather than actual values
  if( length(cross_pred$predvar) != length(temp_quantiles) ){
    stop( paste0("Length of data and cross prediction (with centering) is not equal for ", reg, " region !!!" ) ); 
  }
  
  return( list("MMT"        = MMT, 
               "cross_pred" = cross_pred ) ) 
}


cross_reduced <- function(cb_temp, model, MMT, type="overall", value=NULL){
  "cb_temp: 
   model: 
   MMT:
   type:
   value:
   
   TO DO - ADD OBJECTIVE"
  
  reduced = crossreduce( cb_temp, model,
                         model.link = "log",
                         cen        = MMT,
                         type       = type,
                         value      = value ) ;
  coef = coef( reduced );
  vcov = vcov( reduced );
  
  # omit out
  if( any( !is.finite( coef ) ) ){
    print( "    WARNING: Non-Finite Reduced  Coefficients" );
    print( coef );
  }
  if( any( !is.finite( vcov ) ) ){
    print( "    WARNING: Non-Finite Reduced Covariance Matrix" );
    print( vcov );
  }
  
  return( list("coeff" = coef,
               "vcov"  = vcov  ) )
}


MMT_centered_onebasis <- function(temp_calib, temp_pred, MMT, temp_knots, func_type, 
                                  temp_degree){
  "temp_calib: 
   temp_pred: 
   temp_knots:
   func_type:
   temp_degree:
   
   TO DO - ADD OBJECTIVE"
  
  # temp. basis
  argvar = list( x     = temp_pred, 
                 knots = quantile( temp_calib, temp_knots, na.rm = TRUE ),
                 fun   = func_type,
                 Bound = range( temp_calib, na.rm = TRUE) ) ;
  if( func_type == "bs" ){ argvar =  c(argvar, list( degree = temp_degree )) }
  ob_temp = do.call( onebasis, argvar) ;
  
  # MMT
  argvar_MMT = list( x     = MMT, 
                     knots = quantile( temp_calib, temp_knots, na.rm = TRUE ),
                     fun   = func_type,
                     Bound = range( temp_calib, na.rm = TRUE) ) ;
  if( func_type == "bs" ){ argvar =  c(argvar, list( degree = temp_degree )) }
  ob_temp_MMT = do.call( onebasis, argvar_MMT) ;
  
  # temp. basis centered on MMT
  ob_temp_centered = scale( ob_temp, center = ob_temp_MMT, scale = FALSE );
  return(ob_temp_centered)
}


## from TJ, JB without modifications- 

#  Wald Test of the meta-predictors
##############################################################################
FWALD <- function(model,var){
  ind <- grep(var,names(coef(model)));
  coef <- coef(model)[ind];
  vcov <- vcov(model)[ind,ind];
  waldstat <- coef%*%solve(vcov)%*%coef;
  df <- length(coef);
  return(1-pchisq(waldstat,df));
}


# AIC for overdispersed count
##############################################################################
# THIS FUNCTION RETURNS THE QUASI-AIC STATISTIC FOR GLM MODEL FITTED 
#	THROUGH QUASI-LIKELIHOOD OR THE STANDARD AIC, FOLLOWING THE FORMULA:
#
# QAIC = -2(LOGLIK) + 2*DF*PHI (CFR. PENG 2006 J.R.Stat.Soc.A,162 PAG 190)
#
# THE R OUTPUT DOES NOT PROVIDE BOTH LOGLIK AND AIC FOR QUASI-LIKELIHOOD MODELS
#
# 2 VERSION OF AIC ARE PROVIDED: THE ONE ABOVE AND ONE WHERE THE FIRST TERM 
#	IS REPLACED BY THE DEVIANCE WHEN TYPE="DEV" (COMPARABLE WITH STATA, FOR EXAMPLE)
# THE DIFFERENCE IS ALWAYS AN ADDITIVE CONSTANT (SEE help(extractAIC)), SO THE CHOICE
#	DOESN'T AFFECT THE COMPARISON OF THE MODELS
#
# PHI IS THE OVERDISPERSION PARAMETER, DF THE NUMBER OF NON-ALIASED COEFF 
#
##############################################################################

QAIC <- function(model,type="logLik") {
  QAIC <- vector("list",0)
  
  if(!model$family$family%in%c("poisson","quasipoisson")) {
    stop("only for poisson/quasipoisson family")
  }
  
  phi <- summary(model)$dispersion
  if(type=="dev") {
    QAIC <- deviance(model) + 2*summary(model)$df[3]*phi
  } else {
    loglik <- sum(dpois( model$y, model$fitted.values, log=TRUE))
    QAIC <- -2*loglik + 2*summary(model)$df[3]*phi
  }
  
  return(QAIC)
}

#