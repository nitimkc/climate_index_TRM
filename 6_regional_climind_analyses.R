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
suppressMessages( pacman::p_load(config, profvis, tidyverse, parallel) )

source("climate_index_analyses.R")                 # analyses of climate index in use

CONFIG <- config::get()
print(CONFIG$NAO_THRESHOLD)
print(CONFIG$PARALLEL)


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
n_threshold = n_distinct(NAO$binary_thres, na.rm = TRUE)
thresholds = switch(n_threshold==2, c('pos', 'neg'), NA) # TO DO add what other threshold might be

# threshold_cols = grepl( "thres", colnames(NAO) )
# TO DO -- revise, refactor, automate
# thres_cols_idx = which(grepl( 'thres', colnames(NAO) ))
# dates = NAO$date

# both .rds and parquet are missing date rownames
# TO DO - so many dates variable seem unnecessary, review and refactor
all_dates = as.character( seq(as.Date(CONFIG$NUTS_STARTDATE), as.Date(CONFIG$NUTS_ENDDATE), 1) )
NAO = NAO[NAO$date <= CONFIG$NUTS_ENDDATE, ]

NAO_positive_dates = as.character( NAO[NAO$binary_thres==TRUE, ]$date )
NAO_negative_dates = as.character( NAO[NAO$binary_thres!=TRUE, ]$date )
# TO DO add additional thresholds if any


# ------------------------------------------------------------------------------
# GET ATTRIBUTABLE FRACTION DATA FOR EACH NAO PHASE
# ------------------------------------------------------------------------------

n_simu = CONFIG$N_SIMU
n_groups = length(temp_groups)
n_regions = length(info_region$code)  
col_names = c( sprintf("simu_%s", seq(1:n_simu)), list("attr") )

AF   = array( NA, dim  =    c( n_regions,        n_threshold, 1+n_simu,   n_groups   ),
              dimnames = list( info_region$code, thresholds,  col_names, temp_groups ) )
# ttest_pval = wcox_pval = array( NA, dim  =    c( n_regions,        n_groups   ),
#                                 dimnames = list( info_region$code, temp_groups) )

foldin = paste0( foldout, "AF_ts_simu/" ) 
attr_files = list.files(foldin)
# "AF_ts_simu_pqt/" # parquet files do not have region name
# test = read_parquet(paste0(parqt_folder, parqt_files[200])) 
fname_note = tail( strsplit(sub("\\.[^.]*$", "", CONFIG$NAO_THRESHOLD), "_")[[1]], 1)

start_time = Sys.time() # 10 min parallel/ 24 mins without
if (CONFIG$PARALLEL==TRUE){
  # parallelize
  # ------------
  mean_attr <- function(r, sgn) {
    reg = info_region$code[r]
    print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", reg, ")" ) )
    
    attributions = readRDS(paste0(foldin, reg, ".rds"))
    rownames(attributions) = all_dates
    
    attr_phase = filter(attributions, all_dates %in% as.Date( sgn ) )
    sapply(temp_groups, get_group_mean, attr_phase, simplify="array")
  }
  nbc <- 4
  attr_pos <- mclapply(1:n_regions, mean_attr, NAO_positive_dates, mc.cores=nbc);
  attr_neg <- mclapply(1:n_regions, mean_attr, NAO_negative_dates, mc.cores=nbc);
  for (r in 1:n_regions) {
    AF[r,1,,] <- attr_pos[[r]]
    AF[r,2,,] <- attr_neg[[r]]
  }
  saveRDS(AF, paste0(foldout, "AF_NAO_", fname_note, "_parallelized.rds"))
  
} else {
  # w/out parallelize
  # -----------------
  for (r in 1:n_regions){
    reg = info_region$code[r]
    print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", reg, ")" ) )

    attributions = readRDS(paste0(foldin, reg, ".rds"))
    rownames(attributions) = all_dates
    
    attr_pos_phase = filter(attributions, all_dates %in% NAO_positive_dates)
    attr_neg_phase = filter(attributions, all_dates %in% NAO_negative_dates)
    
    AF[r,1,,] = sapply(temp_groups, get_group_mean, attr_pos_phase, simplify="array")
    AF[r,2,,] = sapply(temp_groups, get_group_mean, attr_neg_phase, simplify="array")
    
    # ttest_pval[r,] = apply( AF[r,,,], 3, function(x) get_pval(x[1,], x[2,] ) )                     # 1x7
    # wcox_pval[r,]  = apply( AF[r,,,], 3, function(x) get_pval(x[1,], x[2,], test_type="Wilcox" ) ) # 1x7
  }
  saveRDS(AF, paste0(foldout, "AF_NAO_", fname_note, ".rds"))
}
end_time = Sys.time()
print(end_time-start_time)


# ------------------------------------------------------------------------------
# COMPUTE STATISTICAL DIFFERENCE OF MEAN BETWEEN EACH NAO PHASE
# ------------------------------------------------------------------------------

ttest_pval = apply( AF, c(1,4), function(x) get_pval(x[1,], x[2,]) )
wcox_pval  = apply( AF, c(1,4), function(x) get_pval(x[1,], x[2,], test_type="Wilcox") )
# TO DO add dimnames

saveRDS(ttest_pval, paste0(foldout, "ttest_pval_AF_NAO_", fname_note, ".rds"))
saveRDS(wcox_pval, paste0(foldout, "wcox_pval_AF_NAO_", fname_note, ".rds"))

print( apply(ttest_pval, 2, FUN=function(x) {length(x[x<=0.05])/length(x)}) )
print( apply(ttest_pval, 2, function(x) table(x<0.05)) )

print( apply(wcox_pval, 2, FUN=function(x) {length(x[x<=0.05])/length(x)}) )
print( apply(wcox_pval, 2, function(x) table(x<0.05)) )


# # ------------------------------------------------------------------------------
# # GROUP BY COUNTRY AND LARGER REGIONS
# # ------------------------------------------------------------------------------
# 
# detrend_type = tail( strsplit(sub("\\.[^.]*$", "", CONFIG$NAO_THRESHOLD), "_")[[1]], 1)
# fname_note   = ifelse(CONFIG$CLIM_IND_SCALED==TRUE, "", detrend_type)
# AF_NAO       = readRDS( paste0(foldout, "AF_NAO", fname_note, ".rds") )
# pval         = readRDS( paste0(foldout, "pval_AF_NAO", fname_note, ".rds") )
# test_sig     = readRDS( paste0(foldout, "test_sig_AF_NAO", fname_note, ".rds") )
# 
# popn = read.table( paste0(CONFIG$ROOT, "indata/", CONFIG$POPULATION), header=TRUE, sep="," )
# popn = popn[popn$year==2024, c(1,3)]
# 
# grps = 1: length(temp_groups)
# simu = 1: 1000
# attr = AF_NAO[,grps,2,simu] - AF_NAO[,grps,1,simu] # simu only
# attr = lapply(seq(dim(attr)[3]), function(x) attr[ , , x])
# 
# for (i in 1:length(attr)){
#   data = as.data.frame(attr[,,i])
#   summarized = popn_wighted_attributions(data)
# }
# 
# plot_df = as.data.frame(attr[,,1]) #;colnames(plot_df) = c("Total")
# 
# AF_pwei_country = popn_wighted_attributions(plot_df, popn, merge_by="location", grp_by="CNTR_CODE")
# 
# test = apply(attr, 3, popn_wighted_attributions)
# 
# 
# popn_wighted_attributions = function(data, popn=popn, merge_by="location", grp_by="CNTR_CODE", cols=temp_groups){
#   data = as.data.frame(data)
#   data = base::merge(popn, data, by.x=merge_by, by.y="row.names")
#   data[grp_by] = substr(data[,merge_by], 1, 2)
# 
#   summarized = data %>%
#     group_by_at(grp_by) %>% # grp_by
#     select(cols) %>%
#     summarise( across(cols,
#                       sum(popn*.x) / sum(popn)
#                       ) ) %>%
#     column_to_rownames(var=grp_by)
#   return(summarized)
# }


