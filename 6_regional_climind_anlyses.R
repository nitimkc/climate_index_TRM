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

# foldname = "AF_ts_simu_pqt/" # parquet is also missing region name
foldname = "AF_ts_simu/"  
foldout_attr = paste0( foldout, foldname ) # TO DO - MOVE FOLDER NAME TO CONFIG
attr_files = list.files(foldout_attr)
# test = read_parquet(paste0(parqt_folder, parqt_files[200])) 

start_time = Sys.time()
for (r in 1:n_regions){
  
  reg = info_region$code[r]
  print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", reg, ")" ) ) 
  attributions = readRDS(paste0(foldout_attr, reg, ".rds"))
  
  rownames(attributions) = all_dates
  
  AF[r,,,] = sapply(temp_groups, get_phases, attributions, NAO_positive_dates, NAO_negative_dates, simplify="array") # 1001x2x7
  ttest_pval[r,] = apply( test, 3, function(x) get_pval(x[,1], x[,2] ) )   # 1x7
  wcox_pval[r,]  = apply( test, 3, function(x) get_pval(x[,1], x[,2], test_type="Wilcox" ) ) # 1x7 
}
end_time = Sys.time()
print(end_time-start_time)

# for (g in 1:n_groups) {
#   
#   g_col = temp_groups[g]
#   row_idx = which(attributions[[g_col]] == TRUE)
#   dates = as.Date( rownames(attributions)[row_idx] )
#   
#   df = attributions[row_idx, 1:(n_simu+1)]
#   if ( nrow(df)>1 ) {
#     rownames(df) = dates
#     positive_df = filter(df, dates %in% NAO_positive_dates)
#     negative_df = filter(df, dates %in% NAO_negative_dates)
#     
#     AF[r,g,1, ] = colMeans(positive_df, na.rm=TRUE)
#     AF[r,g,2, ] = colMeans(negative_df, na.rm=TRUE)
#     
#     # h0 - no difference (equal mean)
#     # reject h0 if pval less than sig level
#     ttest_stat = t.test(AF[r,g,1, ], AF[r,g,2, ], 
#                        var.equal = FALSE,  # diff are norm dist 
#                        paired    = TRUE,   # are each simulations same subject? 
#                        alternative="two.sided") # !=
#     ttest_pval[r,g] = ttest_stat$p.value 
#     
#     wcox_stat = wilcox.test(AF[r,g,1, ], AF[r,g,2, ], 
#                             paired=TRUE,
#                             alternative="two.sided") 
#     wcox_pval[r,g] = ttest_stat$p.value
#   }
# }

# apply(pval, 2, FUN=function(x) {length(x[x<=0.05])/length(x)})
# apply(pval, 2, function(x) table(x<0.05))

scaled = ifelse(CONFIG$CLIM_IND_SCALED==TRUE, "_zscaled", "")
saveRDS(AF, paste0(foldout, "AF_NAO", scaled, "2.rds"))
saveRDS(ttest_pval, paste0(foldout, "ttest_pval_AF_NAO", scaled, ".rds"))
saveRDS(wcox_pval, paste0(foldout, "wcox_pval_AF_NAO", scaled, ".rds"))

# TO DO MOVE WRAPPER FUNCTIONS TO ANOTHER FILE
# ----------------------------------------------

get_phases <- function(temp_grp, data, pos_dates, neg_dates, n_simu=1000){

  row_idx = which(data[[temp_grp]]==TRUE)
  dates = as.Date( rownames(data)[row_idx] )
  data = data[ row_idx, 1:(n_simu+1) ]

  if ( nrow(data)>1 ) {
    rownames(data) = dates
    pos_df = colMeans( filter(data, dates %in% pos_dates), na.rm=TRUE)
    neg_df = colMeans( filter(data, dates %in% neg_dates), na.rm=TRUE)
    combined = cbind( pos_df, neg_df )
  }
  
  return(combined)
}
# test = sapply(temp_groups, get_phases, attributions, NAO_positive_dates, NAO_negative_dates, simplify="array") # 1001x2x7

get_pval <- function(x, y, test_type="Student-t", paired=TRUE, alternative="two.sided"){
  if (test_type=="Student-t"){
    pval = t.test(x, y, var.equal=FALSE, paired=paired, alternative=alternative)$p.value
  } else if (test_type=="Wilcox") {
    pval = wilcox.test(x, y, paired=paired, alternative=alternative)$p.value
  }
  return(pval)
}
# ttest_stat = apply( test, 3, function(x) get_pval(x[,1], x[,2] ) )


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


