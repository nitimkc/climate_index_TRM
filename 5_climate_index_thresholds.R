# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra
# - Epidemiological Model of the Daily EARLY-ADAPT Dataset
# - AUTHOR:
# -   NM as of March 25 2025 
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
suppressMessages( pacman::p_load(config, readxl) )

CONFIG <- config::get()


# ------------------------------------------------------------------------------
# SET PARAMETERS
# ------------------------------------------------------------------------------

# file paths
FOLD_DATA_IN  = paste0( CONFIG$ROOT, "indata/")
FOLD_DATA_CLIM  = paste0(FOLD_DATA_IN, CONFIG$CLIMATE_INDICES)
FOLD_DATA_OUT = paste0( CONFIG$ROOT, "outdata/data_out_daily/")

# attributes
sSEX = CONFIG$SEX
sAGE = CONFIG$AGE
sCAU = CONFIG$CAUSE
iSEN = CONFIG$SENSITIVITY # Change when ready for sensitivity analyses

# attribute specific folder
foldout_name  = paste0(sSEX, "_", sAGE, "_", sCAU, "_iSEN.", iSEN)
foldout       = paste0( FOLD_DATA_OUT, foldout_name, "/" );


# ------------------------------------------------------------------------------
# READ INDEX DATA AND OBTAIN THRESHOLDS
# ------------------------------------------------------------------------------
clim_ind_files = list.files( FOLD_DATA_CLIM )

# print( "Read NAO" )
NAO_idx = which( grepl("NAO", clim_ind_files) )
NAO = read_xlsx( paste0(FOLD_DATA_CLIM, clim_ind_files[NAO_idx]), skip=2 ) 
NAO = NAO |> 
  mutate(date = as.Date(make_datetime(year, month, day)))
  
# # MOV Annex file contains EOF of NAO
# NAO_yearly <- NAO |> 
#   group_by(year) |> 
#   summarize(yearly_avg = mean(nao_index_cdas, na.rm=TRUE)) 
# plot(NAO_yearly$year, NAO_yearly$yearly_avg, type="l", lwd=2)

# NAO[is.na(NAO$nao_index_cdas),]
# NAO = NAO |> 
#   mutate(missing_day = date - lag(date)) |> 
#   filter(missing_day > 1)

# binary threshold where +ve is greater than 0 and -ve otherwise
# ----------------
NAO = NAO |> 
  mutate( scaled_nao = c(scale(nao_index_cdas)) ) |> # -4.017555  3.339121
  mutate( binary_thres = scaled_nao > 0 )
  # mutate( binary_thres = nao_index_cdas>0 )
NAO_thresholds =  NAO[ , c("date", "binary_thres") ]
write.csv( NAO_thresholds, paste0(FOLD_DATA_OUT, "NAO_thresholds.csv"), row.names=FALSE )

# to further filter positive, negative or neutral index. 
# TO DO 
# -----



# ------------------------------------------------------------------------------
# SAVE DETRENDED VERSION
# ------------------------------------------------------------------------------
