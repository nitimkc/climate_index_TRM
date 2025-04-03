# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra (NM)
# - Epidemiological Model of the Daily EARLY-ADAPT Dataset
# - adapted from 
# -   Joan Ballester as of October 2024
# -   Tomás Janos as of February 25 2025
# -   NM as of February 26 2025 
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
suppressMessages( pacman::p_load(config, data.table) )

source("input_mortality_files.R")         # List of Input mortality files
source("preprocessing.R")                 # Check and process read data

CONFIG <- config::get()

select_date <- function(dates, desired_dates, max_lag) {
  
  for (i in 1:length(dates)){
    if (i%%2 != 0) {dates[i] = as.Date(max(dates[i], desired_dates[i]))} # start date
    else           {dates[i] = as.Date(min(dates[i], desired_dates[i]))} # end
  }
  dates[4] = dates[4] - max_lag # no mortality data for last max_lags days 

  return(dates)
}


# ------------------------------------------------------------------------------
# SET PARAMETERS
# ------------------------------------------------------------------------------

# file paths
FOLD_DATA_IN   = paste0( CONFIG$DATAIN )
FOLD_DATA_MAIN = paste0( CONFIG$ROOT, "indata/")
FOLD_DATA_OUT  = paste0( CONFIG$ROOT, "outdata/data_out_daily/")

# attributes
sSEX = CONFIG$SEX
sAGE = CONFIG$AGE
sCAU = CONFIG$CAUSE
iSEN = CONFIG$SENSITIVITY # Change when ready for sensitivity analyses

# attributes specific files and folders
vFILE_MORT     = input_mortality_files( sSEX, sAGE, sCAU )
foldout_name = paste0(sSEX, "_", sAGE, "_", sCAU, "_iSEN.", iSEN)
foldout = paste0( FOLD_DATA_OUT, foldout_name, "/" )
if( !file_test( "-d", foldout ) ){ dir.create( file.path(foldout)); }
# sink( paste0( foldout, "log_", foldout_name, ".txt" ), split = TRUE ) # ?? TO DO - What is sink doing?

# dates
dates_regions = as.Date(c(CONFIG$CALIBRATION_START, CONFIG$CALIBRATION_END,
                          CONFIG$PREDICTION_START,  CONFIG$PREDICTION_END ))
mort_dates = lapply(vFILE_MORT, function(i) select_date(i$DATES, dates_regions, CONFIG$LAG_MAX))
for( i in 1:length(vFILE_MORT) ){ vFILE_MORT[[i]]$DATES = mort_dates[[i]] }

# TO DO - so many dates variable seem unnecessary, review and refactor
dates_NUTS = as.Date(c(CONFIG$CALIBRATION_START, CONFIG$CALIBRATION_END,
                       CONFIG$NUTS_STARTDATE,    CONFIG$NUTS_ENDDATE ))
dates_NUTS[4] = dates_NUTS[4] - CONFIG$LAG_MAX # no mortality data for last max_lags days


# ------------------------------------------------------------------------------
# READ DATA
# ------------------------------------------------------------------------------

start_time = Sys.time() # 40.36237 secs

print( "Read meta table for historical GISCO regions" )
metadata = fread( paste0(FOLD_DATA_IN, CONFIG$META) )

print( "Read temperature" )
temp_hist = fread( paste0( FOLD_DATA_IN, CONFIG$TEMP), )
temp_hist$temp = temp_hist$temp - 273.15; # Kelvin to Celsius

print( paste0("Read mortality for ", length(vFILE_MORT), " regions" ) )
mort_files = lapply(vFILE_MORT, function(i) i$FILENAME)
mort <- rbindlist( lapply(mort_files, fread) )

print( paste0("Read country map") )
country_map = read.table( paste0(FOLD_DATA_MAIN, CONFIG$COUNTRY_MAP), header = TRUE, sep = "," )

end_time = Sys.time()
print(end_time-start_time)


# ------------------------------------------------------------------------------
# PREPROCESS DATA
# ------------------------------------------------------------------------------

start_time = Sys.time() # 8.645938 mins, 25.48462 secs
preprocessed <- data_processing(metadata, temp_hist, mort, foldout, foldout_name,
                                att=list(sex=sSEX, age=sAGE, cause=sCAU), 
                                country_map, vFILE_MORT, CONFIG$LAG_MAX, dates_NUTS)
end_time = Sys.time()
print(end_time-start_time)

calib_mort_temp = preprocessed$calib_mort_temp
all_temp        = preprocessed$all_temp
pred_mort_temp  = preprocessed$pred_mort_temp
pred_temp_MOV   = preprocessed$pred_temp_MOV
info_region     = preprocessed$info_region
info_country    = preprocessed$info_country
calib_NUTS      = preprocessed$calib_NUTS
pred_NUTS       = preprocessed$pred_NUTS
info_NUTS       = preprocessed$info_NUTS
meta_NUTS       = preprocessed$meta_NUTS
meta_regions    = preprocessed$meta_regions

# save
saveRDS( vFILE_MORT,      paste0(foldout, "mort_files_dates") )
saveRDS( dates_NUTS,      paste0(foldout, "dates_NUTS") )

saveRDS( all_temp,        paste0(foldout, "all_temp.rds") ) # temp only for calib and pred
saveRDS( calib_mort_temp, paste0(foldout, "calib_mort_temp.rds") )
saveRDS( pred_mort_temp,  paste0(foldout, "pred_mort_temp.rds") )
saveRDS( pred_temp_MOV,   paste0(foldout, "pred_temp_MOV.rds") )
saveRDS( calib_NUTS,      paste0(foldout, "calib_NUTS.rds") )
saveRDS( pred_NUTS,       paste0(foldout, "pred_NUTS.rds") )

saveRDS( info_region,     paste0(foldout, "info_region.rds") )
saveRDS( info_country,    paste0(foldout, "info_country.rds") )
saveRDS( info_NUTS,       paste0(foldout, "info_NUTS.rds") )
