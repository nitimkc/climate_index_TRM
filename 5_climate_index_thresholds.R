# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra
# - Epidemiological Model of the Daily EARLY-ADAPT Dataset
# - AUTHOR:
# -   NM as of March 25 2025 
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# REQUIRED LIBRARIES, FUNCTIONS AND CONFIG
# ------------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
suppressMessages( pacman::p_load(config, readxl, R.utils, lubridate, splines) )

source("climate_index_analyses.R")                 # analyses of climate index in use

CONFIG <- config::get()


# ------------------------------------------------------------------------------
# SET PARAMETERS
# ------------------------------------------------------------------------------

# file paths
FOLD_DATA_IN    = paste0( CONFIG$ROOT, "indata/")
FOLD_DATA_CLIM  = paste0(FOLD_DATA_IN, CONFIG$CLIMATE_INDICES)
FOLD_DATA_OUT   = paste0( CONFIG$ROOT, "outdata/data_out_daily/")

# attributes
sSEX = CONFIG$SEX
sAGE = CONFIG$AGE
sCAU = CONFIG$CAUSE
iSEN = CONFIG$SENSITIVITY # Change when ready for sensitivity analyses

# attribute specific folder
foldout_name  = paste0( sSEX, "_", sAGE, "_", sCAU, "_iSEN.", iSEN )
foldout       = paste0( FOLD_DATA_OUT, foldout_name, "/" );


# ------------------------------------------------------------------------------
# READ INDEX DATA
# ------------------------------------------------------------------------------
clim_ind_files = list.files( FOLD_DATA_CLIM )

# print( "Read NAO" )
index = "NAO"
file_idx = which( grepl(index, clim_ind_files) )
NAO = read_xlsx( paste0(FOLD_DATA_CLIM, clim_ind_files[file_idx]), skip=2 ) 
n_days = ifelse(leap_year(NAO$date), 366, 365)
NAO = NAO |> 
  mutate( date = as.Date(make_datetime(year, month, day)) ) |> 
  mutate( time = as.numeric(date) ) #|> 
  # mutate( deci_year = as.numeric(year) + (yday(nao$Date)-1) / n_days )

# MOV Annex file contains EOF of NAO
NAO_yearly <- NAO |>
  group_by(year) |>
  summarize(yearly_avg = mean(nao_index_cdas, na.rm=TRUE))
plot(NAO_yearly$year, NAO_yearly$yearly_avg, type="l", lwd=2)

# to further filter positive, negative or neutral index. 
# TO DO 
# -----


# ------------------------------------------------------------------------------
# OBTAIN CLIMATE INDEX THRESHOLD incl from detrended data
# ------------------------------------------------------------------------------

N = length(NAO$nao_index_cdas)
na_rows = is.na(NAO$nao_index_cdas)  
n_decades <- ceiling(length(unique(NAO$year)) / 10)

# detrend_type = CONFIG$DETREND # TO DO ADD to config
detrend_type = "fourier"
if (detrend_type == "none") { 
  # 14126
  # z-scale
  scaled = c(scale(NAO$nao_index_cdas)) 
  NAO$binary_thres = scaled > 0
  
} else if (detrend_type == "ols") { 
  # 14132
  fit_form = nao_index_cdas ~ time
  detrended = get_detrended(NAO, fit_form, "nao_index_cdas", "time" )
  NAO$binary_thres = detrended$binary_threshold
  
} else if (detrend_type == "ns") { 
  # 14142 
  # natural spline of time
  fit_form = nao_index_cdas ~ ns(time, df=n_decades)
  detrended = get_detrended(NAO, fit_form, "nao_index_cdas", "time" )
  NAO$binary_thres = detrended$binary_threshold
  
} else if (detrend_type == "loess") { 
  # 13752 
  fit_form = nao_index_cdas ~ time
  detrended = get_detrended(NAO, fit_form, "nao_index_cdas", "time", loess_span=30) # 30-day local window
  NAO$binary_thres = detrended$binary_threshold
  
} else if (detrend_type == "quadratic") { 
  # 14146
  fit_form = nao_index_cdas ~ poly(time, 2, raw=TRUE)
  detrended = get_detrended(NAO, fit_form, "nao_index_cdas", "time")
  NAO$binary_thres = detrended$binary_threshold
  
} else if (detrend_type == "fourier") { 
  # 13764
  # keep high frequency components only i.e. remove slow trends
  # works best when periodic signals are present.
  x = NAO$year[!na_rows]
  y = NAO$nao_index_cdas[!na_rows]
  interpolate = NAO$year[!na_rows]      # linear interpolation required
  remove_frequency = 0.02               # first 2% of frequencies (adjust as needed)
  
  fourier_filtered = fft( approx( x, y, xout=interpolate )$y ) 
  cutoff = floor(N * remove_frequency)
  cutoff_idx = c(1:cutoff, (N-cutoff+1):N)
  fourier_filtered[cutoff_idx] = 0
  
  detrended = Re( fft(fourier_filtered, inverse=TRUE) / N) # inverse to get high freq signal
  NAO$binary_thres = detrended > 0
  
}
print(sum(NAO$binary_thres, na.rm = TRUE))
write.csv( NAO[ , c("date", "binary_thres") ], 
           paste0(FOLD_DATA_OUT, "NAO_thresholds_", detrend_type, ".csv"), 
           row.names=FALSE )


# ------------------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------------------
ggplot(NAO, aes(x = date)) +
  geom_line(aes(y = nao_index_cdas), color = "blue", alpha = 0.4, linetype = "dashed") +
  geom_line(aes(y = detrended$detrended ), color = "green", alpha=0.5) +
  labs(title = "Detrended NAO",
       y = "NAO Index",
       caption = "Blue: Original | Purple: Fourier Detrended") +
  theme_minimal()





















