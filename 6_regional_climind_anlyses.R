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
suppressMessages( pacman::p_load(config, tidyverse, giscoR, lwgeom, sf, tmap, grid) )

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

NAO = read_csv( paste0(FOLD_DATA_OUT, CONFIG$NAO) )

# TO DO -- revise, refactor, automate

# thres_cols_idx = which(grepl( 'thres', colnames(NAO) ))
# dates = NAO$date

n_threshold = 2
thresholds = c('pos', 'neg')
NAO_positive_dates = filter(NAO, binary_thres==TRUE)[['date']] 
NAO_negative_dates = filter(NAO, binary_thres!=TRUE)[['date']]


# ------------------------------------------------------------------------------
# GET ATTRIBUTATION FRACTION DATA FOR EACH REGION 
# ------------------------------------------------------------------------------

n_simu = 1000
n_groups = length(temp_groups)
n_regions = length(info_region$code)  
col_names = c( sprintf("simu_%s", seq(1:n_simu)), list("attr") )

AF   = array( NA, dim  =    c( n_regions,        n_groups,    n_threshold,  1+n_simu ),
              dimnames = list( info_region$code, temp_groups, thresholds,  col_names ) )
ttest_pval = wcox_pval = array( NA, dim  =    c( n_regions,        n_groups),
                                dimnames = list( info_region$code, temp_groups) )

# both .rds and parquet are missing date rownames
# parquet is also missing the region name
# TO DO - so many dates variable seem unnecessary, review and refactor
all_dates = seq(as.Date(CONFIG$NUTS_STARTDATE), as.Date(CONFIG$NUTS_ENDDATE), 1)

foldout_attr = paste0( foldout, "AF_ts_simu/" ) # TO DO - MOVE FOLDER NAME TO CONFIG
attr_files = list.files(foldout_attr)
# parqt_folder = paste0(foldout, "AF_ts_simu_pqt/")
# parqt_files = list.files( paste0(foldout, "AF_ts_simu_pqt") )
# test = read_parquet(paste0(parqt_folder, parqt_files[200])) 

# ?? Is it faster is sapply is used?

start_time = Sys.time();

for (file in attr_files) {
  
  reg = gsub("\\.rds$", "", file)
  r = which(info_region$code==reg)
  print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", info_region$code[r], ")" ) )
  
  attributions = readRDS( paste0(foldout_attr, file) )
  rownames(attributions) = all_dates
  
  for (g in 1:n_groups) {
    
    g_col = temp_groups[g]
    row_idx = which(attributions[[g_col]] == TRUE)
    dates = as.Date( rownames(attributions)[row_idx] )
    
    df = attributions[row_idx, 1:(n_simu+1)]
    if ( nrow(df)>1 ) {
      rownames(df) = dates
      positive_df = filter(df, dates %in% NAO_positive_dates)
      negative_df = filter(df, dates %in% NAO_negative_dates)
      
      AF[r,g,1, ] = colMeans(positive_df, na.rm=TRUE)
      AF[r,g,2, ] = colMeans(negative_df, na.rm=TRUE)
      
      # h0 - no difference (equal mean)
      # reject h0 if pval less than sig level
      ttest_stat = t.test(AF[r,g,1, ], AF[r,g,2, ], 
                         var.equal = FALSE,  # diff are norm dist 
                         paired    = TRUE,   # are each simulations same subject? 
                         alternative="two.sided") # !=
      ttest_pval[r,g] = ttest_stat$p.value 
      
      wcox_stat = wilcox.test(AF[r,g,1, ], AF[r,g,2, ], 
                              paired=TRUE,
                              alternative="two.sided") 
      wcox_pval[r,g] = ttest_stat$p.value
    }
  }
}

end_time = Sys.time();
print(end_time-start_time);

# apply(pval, 2, FUN=function(x) {length(x[x<=0.05])/length(x)})
# apply(pval, 2, function(x) table(x<0.05))

scaled = ifelse(CONFIG$CLIM_IND_SCALED==TRUE, "_zscaled", "")
saveRDS(AF, paste0(foldout, "AF_NAO", scaled, "2.rds"))
saveRDS(ttest_pval, paste0(foldout, "ttest_pval_AF_NAO", scaled, ".rds"))
saveRDS(wcox_pval, paste0(foldout, "wcox_pval_AF_NAO", scaled, ".rds"))

# ---------------------------------------
# run with sapply ??
# data = attributions
# temp_grp = "Total"
# pos_dates = NAO_positive_dates
# neg_dates = NAO_negative_dates
# 
# get_phases <- function(temp_grp, data, pos_dates, neg_dates, n_simu=1000){
#   
#   row_idx = which(data[[temp_grp]]==TRUE)
#   dates = as.Date( rownames(data)[row_idx] )
#   data = data[ row_idx, 1:(n_simu+1) ]
# 
#   if ( nrow(data)>1 ) {
#     rownames(data) = dates
#     pos_df = colMeans( filter(data, dates %in% pos_dates), na.rm=TRUE)
#     neg_df = colMeans( filter(data, dates %in% neg_dates), na.rm=TRUE)
# 
#     combined = cbind( pos_df,neg_df )
#   }
#   return(combined)
# }
# test = sapply(temp_groups, get_phases, attributions, NAO_positive_dates, NAO_negative_dates, simplify=FALSE) # 7x2x1000
# test_stat = lapply( test, function(x) t.test(x ~ y, var.equal = FALSE, paired = TRUE, alternative="two.sided") )


# ------------------------------------------------------------------------------
# GROUP BY COUNTRY AND LARGER REGIONS
# ------------------------------------------------------------------------------


scaled = ifelse(CONFIG$CLIM_IND_SCALED==TRUE, "_zscaled", "")
AF_NAO = readRDS( paste0(foldout, "AF_NAO", scaled, ".rds") )
pval = readRDS( paste0(foldout, "pval_AF_NAO", scaled, ".rds") )
test_sig = readRDS( paste0(foldout, "test_sig_AF_NAO", scaled, ".rds") )

popn = read.table( paste0(CONFIG$ROOT, "indata/", "eurostat_table_popu_allsex_allage.csv"), 
                   header = TRUE, sep = "," )
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


