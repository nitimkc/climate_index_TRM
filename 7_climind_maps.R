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
suppressMessages( pacman::p_load(config, giscoR, lwgeom, sf, tmap, grid) )

source("maps_EU.R")                 # build maps for specific regions

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
foldout_plot = paste0( FOLD_DATA_OUT, foldout_name, "plots/")
if( !file_test( "-d", foldout ) ){ dir.create( file.path(foldout)); }

# temperature groups
temp_groups = strsplit(CONFIG$TEMP_RANGE, ",")[[1]]

# load required outputs from previous runs
info_region  = readRDS( paste0(foldout, "info_region.rds") )

scaled = ifelse(CONFIG$CLIM_IND_SCALED==TRUE, "_zscaled", "")
AF_NAO = readRDS( paste0(foldout, "AF_NAO", scaled, "2.rds") )
pval = readRDS( paste0(foldout, "ttest_pval_AF_NAO", scaled, ".rds") )
# ttest_pval = readRDS( paste0(foldout, "ttest_pval_AF_NAO", scaled, ".rds") )
sig = ttest_pval <= 0.05
print(colSums(sig, na.rm=TRUE))
colnames(sig) = paste0("sig ", colnames(sig))

# others
# TO DO - ?? MOVE TO CONFIG
# plot_codes = c("MMTV", "TANN","TWIN","TSUM","TP01","TP99","TIQR",
#                "RR99","RR95","RR90","RR75","RR50","RR25","RR10","RR05","RR01",
#                "AFCH","AFTC","AFTH")
# NUTS_plot_codes = plot_codes[! plot_codes %in% c("AFCH","AFTC","AFTH") ]

# bounding box for map to plot
# TO DO - ?? MOVE TO CONFIG
eur_bbox = vec_to_bbox( c(18.5,14.5,59.5,54.0)*10e4 ) # c(24,14.5,74,54)*10e4 // st_bbox( shp_NUTS )
cyprus   = list( "bbox"=c(62.5,15.5,66,18),  "title"="Cyprus" ,            "islands"=c("CY")            )
canaries = list( "bbox"=c(15,9,20.75,12),    "title"="Canaries (Spain)",   "islands"=c("ES701","ES702") ) 
azores   = list( "bbox"=c(9,22,14,28.5),     "title"="Azores (Portugal)",  "islands"=c("PT20")          )
maderia  = list( "bbox"=c(17,14.5,19.75,16), "title"="Maderia (Portugal)", "islands"=c("PT30")          )
other_regions = list(cyprus, canaries, azores, maderia)


# ------------------------------------------------------------------------------
# LOAD MAP OBJECTS AND CREATE BASEMAP
# ------------------------------------------------------------------------------

# copy and cache to tmp dir geojson of all countries in the world
gisco_file = paste0(CONFIG$GISCO, "CNTR_RG_", 
                    CONFIG$RESOLUTION, "M_", 
                    CONFIG$GISCO_YEAR, "_",
                    CONFIG$EPSG_CODE, ".geojson")
file.copy( from = gisco_file, 
           to = gsub( "\\\\", "/", gisco_set_cache_dir(verbose=FALSE) ),
           overwrite = TRUE, recursive = TRUE )
shp_country_all = gisco_get_countries(year=CONFIG$GISCO_YEAR, epsg=CONFIG$EPSG_CODE, 
                                      resolution=CONFIG$RESOLUTION, cache=TRUE, 
                                      update_cache=FALSE, verbose=FALSE) # 257 countries
shp_country = st_read( paste0(CONFIG$DATAIN, CONFIG$COUNTRY_GEO) )       # 51 countries with NUTS info
shp_NUTS    = st_read( paste0(CONFIG$DATAIN, CONFIG$REGION_GEO) )        # regions with NUTS info
shp_NUTS    = get_gisco_codes("region", shp_NUTS, info_region)           # only region within the data

# basemap for all maps to plot
basemap = get_basemap( shp_country_all, eur_bbox, shp_country)
for (or in 1:length(other_regions)){
  or_bbox  = vec_to_bbox( other_regions[[or]]$bbox*10e4 )
  other_regions[[or]]$basemap = get_basemap( shp_country_all, or_bbox, shp_country)
  # print(other_regions[[or]]$basemap)
}

# ------------------------------------------------------------------------------
# MAP ATTRIBUTABLE FRACTION FOR EACH REGION - TOTAL TEMP. GROUP ONLY
# ------------------------------------------------------------------------------

# diff when not scaled, scaled and detrended
start_time = Sys.time()
g = 1
# attr = AF_NAO[,g,,1001]
attr = AF_NAO[,1001,,g]
plot_df = as.data.frame( cbind( attr, "diff"=attr[,2]-attr[,1] ) )
plot_df["sig"] = ifelse(sig[,g]==T, "*", "")
reshp_NUTS = base::merge(shp_NUTS, plot_df, by.x="NUTS_ID", by.y=0)      # add variables to plot

fill_var = "diff"
breaks  = get_breaks(plot_df[,'diff'], plot_code="AFTH")
palette = get_palette(plot_code, breaks) # "brewer.rd_bu" # brewer.blues 
# eur_title = "Attributable Fraction between Positive and Negative phases of NAO"

eur = get_map(basemap, reshp_NUTS, fill_var, breaks, palette, midpoint=NA, alpha=0.7, 
              legend=TRUE, legend_title="", legend_size=0.4, frame=TRUE, sig="sig") +
  tm_shape( shp_country, is.main=FALSE ) +
  tm_borders( lwd = 0.5 )
insets = list()
for (or in 1:length(other_regions)){
  print(other_regions[[or]]$title)
  reshp_NUTS_or = reshp_NUTS[(reshp_NUTS$NUTS_ID %in% other_regions[[or]]$islands), ]
  insets[[or]] = get_map(other_regions[[or]]$basemap, reshp_NUTS_or, fill_var, breaks, 
                         palette, midpoint=NA, alpha=0.7, legend=FALSE, frame=TRUE, 
                         title=other_regions[[or]]$title, sig="sig" )
}
vp_list = list(viewport(x=0.12, y=0.80, width=0.15, height=0.15),
               viewport(x=0.14, y=0.65, width=0.20, height=0.20),
               viewport(x=0.10, y=0.45, width=0.20, height=0.20),
               viewport(x=0.10, y=0.28, width=0.10, height=0.10))
plotfname = paste0(foldout_plot, "AF_NAO", scaled, "_Total_paired.png")
tmap_save(eur, insets_tm=insets, insets_vp=vp_list, filename=plotfname, dpi=600)
end_time = Sys.time()
print(end_time-start_time)


# ------------------------------------------------------------------------------
# MAP ATTRIBUTABLE FRACTION FOR EACH REGION - SIX TEMP. GROUPS
# ------------------------------------------------------------------------------

# all other temp groups
start_time = Sys.time()

grps = 1: length(temp_groups)
# attr = AF_NAO[,grps,2,1001] - AF_NAO[,grps,1,1001] # attr only no simu
attr = AF_NAO[,1001,2,grps] - AF_NAO[,1001,1,grps] # attr only no simu
plot_df = as.data.frame( attr )
plot_df[, colnames(sig)[grps]] = ifelse(sig[,grps]==T, "*", "")
plot_df[, colnames(sig)[grps]][is.na(plot_df[, colnames(sig)[grps]])] = ""
reshp_NUTS = base::merge(shp_NUTS, plot_df, by.x="NUTS_ID", by.y=0)      # add variables to plot

fill_var = c("Total Heat", "Moderate Heat", "Extreme Heat", "Total Cold", "Moderate Cold", "Extreme Cold", "Total")
sig_var = paste0("sig ", fill_var)
breaks  = get_breaks(attr, plot_code="AFTH")
palette = get_palette(plot_code, breaks) # "brewer.rd_bu" # brewer.blues 
# eur_title = "Attributable Fraction between Positive and Negative phases of NAO"

eur = get_map(basemap, reshp_NUTS, fill_var, breaks, palette, midpoint=NA, alpha=0.7, 
              legend=TRUE, legend_title="", legend_size=0.4, frame=TRUE, sig=sig_var)
# insets = list()
# for (or in 1:length(other_regions)){
#   print(other_regions[[or]]$title)
#   reshp_NUTS_or = reshp_NUTS[(reshp_NUTS$NUTS_ID %in% other_regions[[or]]$islands), ]
#   insets[[or]] = get_map(other_regions[[or]]$basemap, reshp_NUTS_or, fill_var, breaks, palette, midpoint=NA,
#                          alpha=0.7, legend=FALSE, frame=TRUE, title=other_regions[[or]]$title )
# }
# vp_list = list(viewport(x=0.12, y=0.80, width=0.15, height=0.15),
#                viewport(x=0.14, y=0.65, width=0.20, height=0.20),
#                viewport(x=0.10, y=0.45, width=0.20, height=0.20),
#                viewport(x=0.10, y=0.28, width=0.10, height=0.10))
plotfname = paste0(foldout_plot, "AF_NAO", scaled, "_Groups_paired.png")
tmap_save(eur, filename=plotfname, height = 8.27, width = 11.69, dpi=600)

end_time = Sys.time()
print(end_time-start_time)


# ------------------------------------------------------------------------------
# POPULATION WEIGHTED AF
# ------------------------------------------------------------------------------

popn = read.table( paste0(CONFIG$ROOT, "indata/", "eurostat_table_popu_allsex_allage.csv"), 
                   header = TRUE, sep = "," )
popn = popn[popn$year==2024, c(1,3)]

grps = 1: length(temp_groups)
n_simu = 1000
attr = AF_NAO[,grps,2,1:n_simu] - AF_NAO[,grps,1,1:n_simu] # simu only

plot_df = as.data.frame(attr[,1,1])
colnames(plot_df) = c("Total")
plot_df = base::merge(popn, plot_df, by.x="location", by.y="row.names")
plot_df["CNTR_CODE"] = substr(plot_df$location, 1, 2)

test = plot_df %>%
  group_by(CNTR_CODE) %>%
  summarise( Total = sum(popu * Total) / sum(popu) ) 

plot_df[, colnames(sig)[grps]] = ifelse(sig[,grps]==T, "*", "")
plot_df[, colnames(sig)[grps]][is.na(plot_df[, colnames(sig)[grps]])] = ""
plot_df["CNTR_CODE"] = substr(rownames(plot_df), 1, 2)
