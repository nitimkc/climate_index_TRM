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

source("climate_index_analyses.R")        # analyses of climate index in use
source("attributable_morality.R")         # attr wrappers

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
conf_int = c(0.025,0.975); # two-tailed 95% Confidence Interval # TO DO move to config


# ------------------------------------------------------------------------------
# READ INDEX DATA AND PREPROCESS
# ------------------------------------------------------------------------------

NAO = read_csv( paste0(FOLD_DATA_OUT, CONFIG$NAO_THRESHOLD) )
n_threshold = n_distinct(NAO$binary_thres, na.rm = TRUE)
phases = c('pos', 'neg')
thresholds = switch(n_threshold==2, phases, NA) # TO DO add what other threshold might be

# TO DO -- revise, refactor, automate
# threshold_cols = grepl( "thres", colnames(NAO) )
# thres_cols_idx = which(grepl( 'thres', colnames(NAO) ))
# dates = NAO$date

# both .rds and parquet are missing date rownames
# TO DO - so many dates variable seem unnecessary, review and refactor
all_dates = as.character( seq(as.Date(CONFIG$NUTS_STARTDATE), as.Date(CONFIG$NUTS_ENDDATE), 1) )
NAO = NAO[NAO$date <= CONFIG$NUTS_ENDDATE, ]

# TO DO add additional thresholds if any
NAO_positive_dates = NAO[NAO$binary_thres==TRUE, ]$date
NAO_pos_djf = as.character( get_period_dates(NAO_positive_dates) )
NAO_pos_ja  = as.character( get_period_dates(NAO_positive_dates, period='summer') )

NAO_negative_dates = NAO[NAO$binary_thres!=TRUE, ]$date
NAO_neg_djf = as.character( get_period_dates(NAO_negative_dates) )
NAO_neg_ja  = as.character( get_period_dates(NAO_negative_dates, period='summer') )

NAO_positive_dates = as.character( NAO[NAO$binary_thres==TRUE, ]$date )
NAO_negative_dates = as.character( NAO[NAO$binary_thres!=TRUE, ]$date )


# ------------------------------------------------------------------------------
# GET ATTRIBUTABLE FRACTION DATA FOR VARIOUS TIME PERIODS
# ------------------------------------------------------------------------------

seasons = c("winter", "summer","all")
n_seasons = length(seasons)
n_simu = CONFIG$N_SIMU
n_groups = length(temp_groups)
n_regions = length(info_region$code)  
col_names = c( sprintf("simu_%s", seq(1:n_simu)), list("attr") )
fname_note = tail( strsplit(sub("\\.[^.]*$", "", CONFIG$NAO_THRESHOLD), "_")[[1]], 1)

AF       = array( NA, dim  =    c( n_regions,        n_threshold, 1+n_simu,     n_groups, n_seasons ),
                  dimnames = list( info_region$code, thresholds,  col_names, temp_groups, seasons ) )
AF_all   = array( NA, dim  =    c( n_regions,         1+n_simu,   n_groups ),
                  dimnames = list( info_region$code, col_names, temp_groups ) )

foldin = paste0( foldout, "AF_ts_simu/" ) 
attr_files = list.files(foldin)
# "AF_ts_simu_pqt/" # parquet files do not have region name
# test = read_parquet(paste0(parqt_folder, parqt_files[200])) 

start = Sys.time() # 14 min parallel/ 24 mins without

if (CONFIG$PARALLEL==TRUE){
  # parallelize
  # ------------
  nbc <- 8
  mean_attr <- function(r, sgn) {
    reg = info_region$code[r]
    # print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", reg, ")" ) )
    
    attributions = readRDS(paste0(foldin, reg, ".rds"))
    rownames(attributions) = all_dates
    attr_phase = filter(attributions, all_dates %in% sgn )
    sapply(temp_groups, get_group_mean, attr_phase, simplify="array")
  }
  
  # seasons of each phases
  start_time = Sys.time()
  attr_pos_djf <- mclapply(1:n_regions, mean_attr, NAO_pos_djf, mc.cores=nbc);
  attr_neg_djf <- mclapply(1:n_regions, mean_attr, NAO_neg_djf, mc.cores=nbc);
  print("Calculation for winter complete")
  end_time = Sys.time()
  print(end_time-start_time)
  
  start_time = Sys.time()
  attr_pos_ja  <- mclapply(1:n_regions, mean_attr, NAO_pos_ja, mc.cores=nbc);
  attr_neg_ja  <- mclapply(1:n_regions, mean_attr, NAO_neg_ja, mc.cores=nbc);
  print("Calculation for summer complete")
  end_time = Sys.time()
  print(end_time-start_time)
  
  # each phases
  start_time = Sys.time()
  attr_pos     <- mclapply(1:n_regions, mean_attr, NAO_positive_dates, mc.cores=nbc);
  attr_neg     <- mclapply(1:n_regions, mean_attr, NAO_negative_dates, mc.cores=nbc);
  print("Calculation for each phase complete")
  end_time = Sys.time()
  print(end_time-start_time)
  
  # whole period (both phases)
  start_time = Sys.time()
  attr_all <- mclapply(1:n_regions, mean_attr, all_dates, mc.cores=nbc);
  print("Calculation for entire period complete")
  end_time = Sys.time()
  print(end_time-start_time) 
  
  for (r in 1:n_regions) {
    # whole period
    AF_all[r,,] <- attr_all[[r]]
    
    # positive
    AF[r,1,,,1] <- attr_pos_djf[[r]]
    AF[r,1,,,2] <- attr_pos_ja[[r]]
    AF[r,1,,,3] <- attr_pos[[r]]
    
    # negative
    AF[r,2,,,1] <- attr_neg_djf[[r]]
    AF[r,2,,,2] <- attr_neg_ja[[r]]
    AF[r,2,,,3] <- attr_neg[[r]]
  }
  saveRDS(AF_all, paste0(foldout, "AF_NAO_", fname_note, "_wholeperiod", "_parallelized.rds"))
  saveRDS(AF,     paste0(foldout, "AF_NAO_", fname_note, "_seasons", "_parallelized.rds"))
  
} else {
  # w/out parallelize
  # -----------------
  # test r = 1
  for (r in 1:n_regions) {
    reg = info_region$code[r]
    print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", reg, ")" ) )
    
    attributions = readRDS(paste0(foldin, reg, ".rds"))
    rownames(attributions) = all_dates
    
    # seasons of each phase
    attr_pos_djf   = filter(attributions, all_dates %in% NAO_pos_djf)
    attr_neg_djf   = filter(attributions, all_dates %in% NAO_neg_djf)
    attr_pos_ja    = filter(attributions, all_dates %in% NAO_pos_ja)
    attr_neg_ja    = filter(attributions, all_dates %in% NAO_neg_ja)
    
    # each phase
    attr_pos = filter(attributions, all_dates %in% NAO_positive_dates)
    attr_neg = filter(attributions, all_dates %in% NAO_negative_dates)
    
    # whole period (both phases) 
    start_time = Sys.time()
    AF_all[r,,] = sapply(temp_groups, get_group_mean, attributions, simplify="array")
    print("Calculation for entire period complete")
    end_time = Sys.time()
    print(end_time-start_time)
    
    start_time = Sys.time()
    # positive
    AF[r,1,,,1] = sapply(temp_groups, get_group_mean, attr_pos_djf, simplify="array")
    AF[r,1,,,2] = sapply(temp_groups, get_group_mean, attr_pos_ja, simplify="array")
    AF[r,1,,,3] = sapply(temp_groups, get_group_mean, attr_pos, simplify="array")
    print("Calculation for positive phases complete")
    end_time = Sys.time()
    print(end_time-start_time)
    
    start_time = Sys.time()
    # negative
    AF[r,2,,,1] = sapply(temp_groups, get_group_mean, attr_neg_djf, simplify="array")
    AF[r,2,,,2] = sapply(temp_groups, get_group_mean, attr_neg_ja, simplify="array")
    AF[r,2,,,3] = sapply(temp_groups, get_group_mean, attr_neg, simplify="array")
    print("Calculation for negative phases complete")
    end_time = Sys.time()
    print(end_time-start_time)
  }
  saveRDS(AF_all, paste0(foldout, "AF_NAO_", fname_note, "_wholeperiod", ".rds"))
  saveRDS(AF,     paste0(foldout, "AF_NAO_", fname_note, "_seasons", ".rds"))
}
end = Sys.time()
print(end-start)

# summaries and confidence intervals by regions and countries
# -----------------------------------------------------------

# AF_all = readRDS(paste0(foldout, "AF_NAO_", fname_note, "_wholeperiod_parallelized", ".rds"))
# AF = readRDS(paste0(foldout, "AF_NAO_", fname_note, "_seasons_parallelized", ".rds"))
  
# for whole prediction period
AF_all = apply(AF_all, 3, get_geo_groups, info_region, simplify=FALSE)
AF_confint_all = sapply(AF_all, get_confidence_interval, simplify=FALSE)
AF_confint_all = simplify2array(AF_confint_all)
saveRDS(AF_confint_all, paste0(foldout, "AF_confint_predperiod", ".rds"))
print("confidence intervals at group-level complete")

# for NAO phases and seasons within the prediction period
phase_list = list()
for (phs in (1:n_threshold)){
  
  season_list = list()
  for (ssn in (1:n_seasons)) {
    season = seasons[ssn]
    attr = apply(AF[,phs,,,ssn], 3, get_geo_groups, info_region, simplify=FALSE)
    attr_confint = sapply(attr, get_confidence_interval, simplify=FALSE)
    attr_confint = simplify2array(attr_confint)
    
    season_list = append(season_list, list(attr_confint))
  }
  AF_confint_seasons = do.call(abind::abind, append(season_list, list(rev.along = 0)))
  dimnames(AF_confint_seasons)[[length(dimnames(AF_confint_seasons))]] = seasons
  
  phase_list = append(phase_list, list(AF_confint_seasons))
}
AF_confint_seasons_phases = do.call(abind::abind, append(phase_list, list(rev.along = 0)))
dimnames(AF_confint_seasons_phases)[[length(dimnames(AF_confint_seasons_phases))]] = phases
saveRDS(AF_confint_seasons_phases, paste0(foldout, "AF_NAO_", fname_note, "_confint_phases_seasons.rds"))
print("confidence intervals for phases and seasons complete")

# AF2_confint_seasons_phases = readRDS(paste0(foldout, "AF_NAO_", fname_note, "_confint_phases_seasons.rds"))

# statistical significance test
# -----------------------------
ttest_pval = apply( AF, c(1,4,5), function(x) get_pval(x[1,], x[2,]) )
wcox_pval  = apply( AF, c(1,4,5), function(x) get_pval(x[1,], x[2,], test_type="Wilcox") )
# TO DO add dimnames
print("statistical tests complete")

saveRDS(ttest_pval, paste0(foldout, "ttest_pval_AF_NAO_", fname_note, "_seasons.rds"))
saveRDS(wcox_pval, paste0(foldout, "wcox_pval_AF_NAO_", fname_note, "_seasons.rds"))

print( apply(ttest_pval, c(2,3), FUN=function(x) {length(x[x<=0.05])/length(x)}) )
print( apply(ttest_pval, c(2,3), function(x) sum(x<0.05)) )

print( apply(wcox_pval, c(2,3), FUN=function(x) {length(x[x<=0.05])/length(x)}) )
print( apply(wcox_pval, c(2,3), function(x) sum(x<0.05)) )


# ------------------------------------------------------------------------------
# POPULATION WEIGHTED AF/AN
# ------------------------------------------------------------------------------

# grp = "Total"
# data = AF_all[,,grp]
# dim(data)
# test = plot_df %>%
#   group_by(CNTR_CODE) %>%
#   summarise( Total = sum(popu * Total) / sum(popu) ) 


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


