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
suppressMessages( pacman::p_load(config, tidyverse, giscoR, lwgeom, sf, tmap, grid) )

source("climate_index_analyses.R")                 # analyses of climate index in use

CONFIG <- config::get()


# ------------------------------------------------------------------------------
# SET PARAMETERS
# ------------------------------------------------------------------------------

# file paths
FOLD_DATA_IN     = paste0( CONFIG$ROOT, "indata/")
FOLD_DATA_OUT    = paste0( CONFIG$ROOT, "outdata/data_out_daily/")

# attributes
sSEX = CONFIG$SEX
sAGE = CONFIG$AGE
sCAU = CONFIG$CAUSE
iSEN = CONFIG$SENSITIVITY # Change when ready for sensitivity analyses

# attribute specific folder
foldout_name = paste0(sSEX, "_", sAGE, "_", sCAU, "_iSEN.", iSEN, "/")
foldout = paste0( FOLD_DATA_OUT, foldout_name)

# load required outputs from previous runs
info_region  = readRDS( paste0(foldout, "info_region.rds") )

# temperature groups
temp_groups = strsplit(CONFIG$TEMP_RANGE, ",")[[1]]

# ------------------------------------------------------------------------------
# READ INDEX DATA AND PREPROCESS
# ------------------------------------------------------------------------------

NAO = read_csv( paste0(FOLD_DATA_OUT, CONFIG$NAO_THRESHOLD) )
n_threshold = 2
thresholds = c('pos', 'neg')
NAO_positive_dates = filter(NAO, binary_thres==TRUE)[['date']] 
NAO_negative_dates = filter(NAO, binary_thres!=TRUE)[['date']]

threshold_cols = grepl( "thres", colnames(NAO) )
# TO DO -- revise, refactor, automate
# thres_cols_idx = which(grepl( 'thres', colnames(NAO) ))
# dates = NAO$date


# ------------------------------------------------------------------------------
# GET ATTRIBUTABLE FRACTION DATA FOR EACH NAO PHASE
# ------------------------------------------------------------------------------

n_simu = 1000
n_groups = length(temp_groups)
n_regions = length(info_region$code)  
col_names = c( sprintf("simu_%s", seq(1:n_simu)), list("attr") )

AF   = array( NA, dim  =    c( n_regions,        1+n_simu,  n_threshold, n_groups    ),
              dimnames = list( info_region$code, col_names, thresholds,  temp_groups ) )
ttest_pval = wcox_pval = array( NA, dim  =    c( n_regions,        n_groups),
                                dimnames = list( info_region$code, temp_groups) )

# both .rds and parquet are missing date rownames
# TO DO - so many dates variable seem unnecessary, review and refactor
all_dates = seq(as.Date(CONFIG$NUTS_STARTDATE), as.Date(CONFIG$NUTS_ENDDATE), 1)

foldin = paste0( foldout, "AF_ts_simu/" ) 
attr_files = list.files(foldin)
# "AF_ts_simu_pqt/" # parquet files do not have region name
# test = read_parquet(paste0(parqt_folder, parqt_files[200])) 

start_time = Sys.time() # 1.038734 hours
for (r in 1:n_regions){
  
  reg = info_region$code[r]
  print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", reg, ")" ) ) 
  attributions = readRDS(paste0(foldin, reg, ".rds"))
  rownames(attributions) = all_dates
  
  AF[r,,,] = sapply(temp_groups, get_phases, attributions, 
                    NAO_positive_dates, NAO_negative_dates, 
                    simplify="array") # 1001x2x7
  ttest_pval[r,] = apply( AF[r,,,], 3, function(x) get_pval(x[,1], x[,2] ) )                     # 1x7
  wcox_pval[r,]  = apply( AF[r,,,], 3, function(x) get_pval(x[,1], x[,2], test_type="Wilcox" ) ) # 1x7 
}
end_time = Sys.time()
print(end_time-start_time)

print( apply(ttest_pval, 2, FUN=function(x) {length(x[x<=0.05])/length(x)}) )
print( apply(ttest_pval, 2, function(x) table(x<0.05)) )

print( apply(wcox_pval, 2, FUN=function(x) {length(x[x<=0.05])/length(x)}) )
print( apply(wcox_pval, 2, function(x) table(x<0.05)) )

detrend_type = tail( strsplit(sub("\\.[^.]*$", "", CONFIG$NAO_THRESHOLD), "_")[[1]], 1)
fname_note = ifelse(CONFIG$CLIM_IND_SCALED==TRUE, "", detrend_type)
saveRDS(AF, paste0(foldout, "AF_NAO", fname_note, ".rds"))
saveRDS(ttest_pval, paste0(foldout, "ttest_pval_AF_NAO", fname_note, ".rds"))
saveRDS(wcox_pval, paste0(foldout, "wcox_pval_AF_NAO", fname_note, ".rds"))


# ------------------------------------------------------------------------------
# GROUP BY COUNTRY AND LARGER REGIONS
# ------------------------------------------------------------------------------

detrend_type = tail( strsplit(sub("\\.[^.]*$", "", CONFIG$NAO_THRESHOLD), "_")[[1]], 1)
fname_note   = ifelse(CONFIG$CLIM_IND_SCALED==TRUE, "", detrend_type)
AF_NAO       = readRDS( paste0(foldout, "AF_NAO", fname_note, ".rds") )
pval         = readRDS( paste0(foldout, "pval_AF_NAO", fname_note, ".rds") )
test_sig     = readRDS( paste0(foldout, "test_sig_AF_NAO", fname_note, ".rds") )

popn = read.table( paste0(CONFIG$ROOT, "indata/", CONFIG$POPULATION), header=TRUE, sep="," )
popn = popn[popn$year==2024, c(1,3)]

grps = 1: length(temp_groups)
simu = 1: 1000
attr = AF_NAO[,grps,2,simu] - AF_NAO[,grps,1,simu] # simu only
attr = lapply(seq(dim(attr)[3]), function(x) attr[ , , x])

for (i in 1:length(attr)){
  data = as.data.frame(attr[,,i])
  summarized = popn_wighted_attributions(data)
}

plot_df = as.data.frame(attr[,,1]) #;colnames(plot_df) = c("Total")

AF_pwei_country = popn_wighted_attributions(plot_df, popn, merge_by="location", grp_by="CNTR_CODE")

test = apply(attr, 3, popn_wighted_attributions)


popn_wighted_attributions = function(data, popn=popn, merge_by="location", grp_by="CNTR_CODE", cols=temp_groups){
  data = as.data.frame(data)
  data = base::merge(popn, data, by.x=merge_by, by.y="row.names")
  data[grp_by] = substr(data[,merge_by], 1, 2)

  summarized = data %>%
    group_by_at(grp_by) %>% # grp_by
    select(cols) %>%
    summarise( across(cols, 
                      sum(popn*.x) / sum(popn) 
                      ) ) %>% 
    column_to_rownames(var=grp_by)
  return(summarized)
}


