# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra
# - Check and preprocess temperature, mortality and meta tables
# - adapted from 
# -   Joan Ballester as of October 2024
# -   Tomás Janos as of February 25 2025
# -   NM as of February 27 2025 
# ------------------------------------------------------------------------------


if (!require("pacman")) install.packages("pacman")
suppressMessages( pacman::p_load(dplyr, lubridate) )


missing_items <- function(v1, v2){
  "v1 and v2 are vector of at least one non NA elements.
   
   Check for items in v1 not existing in v2"
  
  missing = v1 %in% v2
  if( any(!missing) ){ stop(paste0("missing in v2 - ", list(v1[!missing]) )) }
}


filter_timeperiod <- function(data, dates, max_lag){
  "data: A table of data to filter with date column to satisfy filter conditions
   dates: A vector of maximum and minimum date for filter conditions
   max_lag: An integer to add to the maximum date (used mainly for prediction period)
  
   Filter 'data' to desired time period based on min and max dates in 'dates' argument."
  
  data = data[ which( data$date >= dates[1] & data$date <= dates[2] + max(max_lag,0) ), ]
  
  if( is.finite( dates[1] ) & any( data$date < dates[1] ) )
  { stop(paste0( "Start date out of range!!!" ) )}
  if( is.finite( dates[2] ) & any( data$date > dates[2] + max(max_lag,0) ) )
  { stop(paste0( "Final date out of range!!!" ) )}
  if( any( diff( data$date ) != 1 ) )
  { stop( paste0( "Missing dates in calibration period" ) )}
  
  return (data)
}


data_processing <- function( metadata, temp_hist, mort, savepath, savefilename,
                             att, country_map, vFILE_MORT, max_lag, dates_NUTS, popn){
  "Data processing is specific to the data files in use
  1. Metadata 
  2. Historical temperature
  3. Mortality
  4. Country map
  Refer to config.yml for input data details."
  
  start_time = Sys.time()
  
  print ("Initial data check.")
  # ---------------------------
  
  # all regions included in each data table
  regions_meta = metadata$location
  regions_unq_temp = unique( temp_hist$location )
  regions_unq_mort = unique( mort     $location )
  
  missing_items(regions_meta, regions_unq_temp)
  missing_items(regions_unq_temp, regions_meta)
  missing_items(regions_unq_mort, regions_unq_temp)
  # we have temperature for regions without mortality data 
  
  # sex, age and cause as intended
  if( any( mort$sex   != att$sex ) )  { stop(paste0("Invalid sex attribute in mortality table - ", unique(mort$sex) ) )}
  if( any( mort$age   != att$age ) )  { stop(paste0("Invalid age attribute in mortality table - ", unique(mort$age) ) )}
  if( any( mort$cause != att$cause ) ){ stop(paste0("Invalid mortality cause attribute in mortality table - ", unique(mort$cause) ) )}
  mort$sex = mort$age = mort$cause = NULL; # ?? TO DO - Check why is this set to NULL instead of being discarded
  
  # any non-finite temp. & mort. (negative) values
  finite_temp = is.finite( temp_hist$temp )
  if( any( !finite_temp ) ) { 
    print( temp_hist[ which( !finite_temp ), ] )
    stop( paste0(sum( !finite_temp ), " non-finite temperature values found.") )
  }
  
  finite_mort = is.finite( mort$mort )
  if( any( !finite_mort ) ) {
    print( mort[ which( !finite_mort ), ] )
    print( paste0("WARNING for Non-Finite Mortality Count(s): ", sum( !finite_mort )) )
  }
  if( any( finite_mort & mort$mort < 0 ) ) {
    print( mort[ which( finite_mort & mort$mort < 0 ), ] );
    stop( paste0("Negative Mortality Counts: ", sum(finite_mort & mort$mort < 0)) )
  }
  
  # any out-of-range temp. values
  min_temp = -50; max_temp = +50
  if( any( min_temp > temp_hist$temp | temp_hist$temp > max_temp ) ) {
    print( temp_hist[ which( min_temp > temp_hist$temp | temp_hist$temp > max_temp ), ] )
    stop( paste0(sum( min_temp > temp_hist$temp | temp_hist$temp > max_temp ), " out-of-range temperature values found.") )
  }
  
  print ("Merge temperature and mortality data.")
  # ---------------------------------------------
  
  # obtain temp. for regions with mort. data
  mort_temp     = base::merge(mort, temp_hist, by=c("location","date"), all.x=TRUE)
  mort_temp     = mort_temp[ order(mort_temp$location, mort_temp$date), ]
  mort_temp$dow = lubridate::wday( mort_temp$date, week_start = 1 )   # add day of the week
  
  # check all regions in merged data are included in meta table
  metadata_unq_regions = unique(metadata$location)
  mort_temp_unq_regions = unique(mort_temp$location)
  missing_items( mort_temp_unq_regions, metadata_unq_regions )
  
  # keep meta data of regions in merged data only
  metadata = metadata[ metadata$location %in% mort_temp_unq_regions, ]
  
  # save( metadata, mort_temp, file = paste0( savepath, "ORIGINAL_DATATABLES_", savefilename, ".RData" ) )
  
  # # TO DO - If needed understand why the following is required
  # # remove regions for which the fit does not converge in certain cases
  # # ADD ARGUMENTS TO THE FUNCTION TO MAKE THIS WORK LOCALLY
  # # LATER TODO move this with data_processing function
  # if( sSEX == "allsex" & sAGE == "age085ppp" & sCAU == "allcause" ){
  #   nonCONVERGED_age085ppp = strsplit(NONCONVERGED_age085ppp, ",")[[1]]
  #   DATATABLE_DATA = DATATABLE_DATA[ !DATATABLE_DATA$location %in% nonCONVERGED_age085ppp, ]
  #   TABLE_INFO_DATA = TABLE_INFO_DATA[ !TABLE_INFO_DATA$location %in% nonCONVERGED_age085ppp, ]
  # }
  # if( sSEX == "allsex" & sAGE == "allage" & sCAU == "j00j99" ){
  #   nonCONVERGED_j00j99 <- strsplit(NONCONVERGED_j00j99, ",")[[1]]
  #   DATATABLE_DATA = DATATABLE_DATA[ !DATATABLE_DATA$location %in% nonCONVERGED_j00j99, ]
  #   TABLE_INFO_DATA = TABLE_INFO_DATA[ !TABLE_INFO_DATA$location %in% nonCONVERGED_j00j99, ]
  # }
  
  print ("Regions and countries with mortality - ")
  # -----------------------------------------------
  
  # meta data of regions to be called info_regions
  region = list("code" = metadata$location, #unique
                "name" = dplyr::if_else( is.na(metadata$location_eng),
                                         metadata$location_off,
                                         metadata$location_eng ),
                "country_code" = metadata$nuts0)
  region$eea_subregion <- country_map$EEA_SUBREGION[ match(region$country_code, country_map$CODE) ]
  region$popn24 <- popn$popu[ match(popn$location, region$code) ]
  
  # meta data of countries to be called info_country
  country_code = unique(metadata$nuts0)
  country = list("code" = country_code,
                 "name" = country_map$NAME[ match(country_code, country_map$CODE) ],
                 "eea_subregion" = country_map$EEA_SUBREGION[ match(country_code, country_map$CODE) ]
                 )
  
  print( paste0("        Number of deaths: ", sum(mort_temp$mort, na.rm=TRUE)) )
  print( paste0("        Number of regions: ", length(region$code)) )
  print( paste0("        Number of EEA sub-regions: ", length(unique(region$eea_subregion))) )
  print( paste0("        Number of countries: ", length(country$code)) )
  
  # obtain calibration data
  #   1. by region
  calib_mort_temp = split(mort_temp, by="location", keep.by=FALSE)
  calib_mort_temp = pred_mort_temp = calib_mort_temp[region$code]
  
  for( r in region$code ){
    pattern = substr(r, 1,3)
    if     ( pattern %in% paste0( "UK", LETTERS[ 3:12] ) ){ file_idx = which( names(vFILE_MORT) == "UKEW" ); }
    else if( pattern %in% paste0( "UK", LETTERS[ 3:11] ) ){ file_idx = which( names(vFILE_MORT) == "UKE"  ); } 
    else if( pattern %in% paste0( "UK", LETTERS[12:12] ) ){ file_idx = which( names(vFILE_MORT) == "UKW"  ); }
    else if( pattern %in% paste0( "UK", LETTERS[13:13] ) ){ file_idx = which( names(vFILE_MORT) == "UKS"  ); }
    else if( pattern %in% paste0( "UK", LETTERS[14:14] ) ){ file_idx = which( names(vFILE_MORT) == "UKNI" ); }
    else { file_idx = which( names(vFILE_MORT) == substr(r, 1,2) ) }
    if( length(file_idx) != 1 ) { stop( paste0("Invalid file index ", file_idx, " in Input Mortality File.") ); }
    
    #   2. by period
    dates_mort = lapply( vFILE_MORT, function(i) i$DATES )  # move to main script once use of vFILE_MORT above is determined
    calib_mort_temp[[r]] = filter_timeperiod( calib_mort_temp[[r]], dates_mort[[file_idx]][1:2], max_lag )
    pred_mort_temp[[r]]  = filter_timeperiod( pred_mort_temp[[r]],  dates_mort[[file_idx]][3:4], max_lag )
  }
  
  # An exhaustive ordered vector with all possible temperature values for MMT calculation
  #   done to minimize the occurrence of negative predicted attributable fractions 
  #   due to lack of MMT precision
  # // NM obtain association for each existing temperature value for a given region for precise
  # // MMT calculation and not miss some temperature values for which the MMT in turn might exist
  all_temp = mapply(function(x, y) sort(unique(c(x$temp, y$temp))), calib_mort_temp, pred_mort_temp)
  
  # obtain prediction data
  pred_temp_MOV = temp_hist[ temp_hist$location %in% region$code, ] # 1. by region 
  pred_temp_MOV = split(pred_temp_MOV, by="location", keep.by=FALSE)
  pred_temp_MOV = pred_temp_MOV[region$code]
  pred_temp_MOV = lapply(pred_temp_MOV, filter_timeperiod, dates_NUTS[3:4], max_lag=21) # 2. by period (temperature without mortality)
  
  print ("NUTS 0-3 regions - ")
  # ---------------------------
  
  # meta data of NUTS 0-3 regions in the last GISCO year to be called info_NUTS
  metadata_NUTS = arrange( metadata[which(metadata$year_shapefile == max(metadata$year_shapefile)
                                          & metadata$level >= 0), ],
                           level )
  NUTS = list("code" = metadata_NUTS$location,
              "name" = dplyr::if_else( is.na( metadata_NUTS$location_eng ),
                                       metadata_NUTS$location_off,
                                       metadata_NUTS$location_eng ))
  print( paste0("        Number of NUTS regions: ", length(NUTS$code)) )
  
  # obtain NUTS calibration and prediction data
  #   1. by region
  temp_NUTS = temp_hist[ temp_hist$location %in% NUTS$code, ]
  
  calib_NUTS = pred_NUTS = split(temp_NUTS, by="location", keep.by=FALSE)
  names(calib_NUTS) = names(pred_NUTS) = NUTS$code
  
  #   2. by period
  calib_NUTS = lapply(calib_NUTS, filter_timeperiod, dates_NUTS[1:2], max_lag=21)
  pred_NUTS  = lapply(pred_NUTS,  filter_timeperiod, dates_NUTS[3:4], max_lag=21)
  
  
  return( list("calib_mort_temp" = calib_mort_temp, 
               "all_temp"        = all_temp,
               "pred_mort_temp"  = pred_mort_temp,
               "pred_temp_MOV"   = pred_temp_MOV,
               "calib_NUTS"      = calib_NUTS,
               "pred_NUTS"       = pred_NUTS, 
               "meta_regions"    = metadata, 
               "meta_NUTS"       = metadata_NUTS, 
               "info_region"     = region, 
               "info_country"    = country, 
               "info_NUTS"       = NUTS ) )
}












