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
suppressMessages( pacman::p_load(config, giscoR, lwgeom, sf, tmap) )

source("NUTS_EU_maps.R")                 # Check and process read data

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
info_country = readRDS( paste0(foldout, "info_country.rds") )
info_NUTS = readRDS( paste0(foldout, "info_NUTS.rds") )

# others
temp_groups = strsplit(CONFIG$TEMP_RANGE, ",")[[1]]


# ------------------------------------------------------------------------------
# READ INDEX DATA AND PREPROCESS
# ------------------------------------------------------------------------------

NAO = read_csv( paste0(FOLD_DATA_OUT, CONFIG$NAO) )
# thres_cols_idx = which(grepl( 'thres', colnames(NAO) ))
# dates = NAO$date
# TO DO -- revise, refactor, automate

n_threshold = 2
thresholds = c('pos', 'neg')
NAO_positive_dates = filter(NAO, binary_thres==TRUE)[['date']] 
NAO_negative_dates = filter(NAO, binary_thres!=TRUE)[['date']]


# ------------------------------------------------------------------------------
# READ ATTRIBUTATION FRACTION DATA FOR EACH REGION 
# ------------------------------------------------------------------------------

n_simu = 1000
n_groups = length(temp_groups)
n_regions = length(info_region$code)  
col_names = c( sprintf("simu_%s", seq(1:n_simu)), list("attr") )

AF     = array( NA, dim  =    c( n_regions,        n_groups,    n_threshold,  1+n_simu ),
                dimnames = list( info_region$code, temp_groups, thresholds,  col_names ) )
t_stat = array( NA, dim  =    c( n_regions,        n_groups),
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

for (file in attr_files[299:654]) {
  
  reg = gsub("\\.rds$", "", file)
  r = which(info_region$code==reg)
  print( paste0( "  Region ", r, " / ", n_regions, ": ", info_region$name[r], " (", info_region$code[r], ")" ) )
  
  attributions = readRDS( paste0(foldout_attr, file) )
  rownames(attributions) = all_dates
  
  for (g in 1:n_groups) {
    g_col = temp_groups[g]
    
    df_idx = which(attributions[[g_col]] == TRUE)
    dates = as.Date( rownames(attributions)[df_idx] )
    
    df = attributions[df_idx, 1:(n_simu+1)]
    if ( nrow(df)>1 ) {
      rownames(df) = dates
      
      positive_df = filter(df, dates %in% NAO_positive_dates)
      negative_df = filter(df, dates %in% NAO_negative_dates)
      
      AF[r,g,1, ] = colMeans(positive_df, na.rm=TRUE); AF[r,g,2, ] = colMeans(negative_df, na.rm=TRUE)
      test_stat = t.test(AF[r,g,1, ], AF[r,g,2, ], var.equal=TRUE, alternative="two.sided")
      t_stat[r,g] = test_stat$p.value
      print(paste0(temp_groups[g], " ",  test_stat$p.value > 0.05))
    }
  }
}

saveRDS(AF, paste0(foldout, "AF_NAO.rds"))
saveRDS(t_stat, paste0(foldout, "t_stat_AF_NAO.rds"))
# AF_NAO = readRDS( paste0(foldout, "AF_NAO.rds") )
# t_stat = readRDS( paste0(foldout, "t_stat_AF_NAO.rds") )


# ------------------------------------------------------------------------------
# LOAD MAP OBJECTS
# ------------------------------------------------------------------------------

plot_codes = c("MMTV", "TANN","TWIN","TSUM","TP01","TP99","TIQR",
               "RR99","RR95","RR90","RR75","RR50","RR25","RR10","RR05","RR01",
               "AFCH","AFTC","AFTH")
NUTS_plot_codes = plot_codes[! plot_codes %in% c("AFCH","AFTC","AFTH") ]

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
                                      update_cache=FALSE, verbose=FALSE)
shp_country_EU = st_read( paste0(CONFIG$DATAIN, CONFIG$COUNTRY_GEO) )

shp_NUTS_EU    = st_read( paste0(CONFIG$DATAIN, CONFIG$REGION_GEO) )
shp_NUTS_EU    = get_gisco_codes("region", shp_NUTS_EU, info_region)
# gisco_country = get_gisco_codes("country", shp_NUTS_EU, info_country)

plot_code = "RR99"
plot_df = AF_NAO[ , 1, 1, 1001] # Total, positive and attr
reshp_NUTS_EU = base::merge(shp_NUTS_EU, plot_df, by.x="NUTS_ID", by.y=0)
breaks  = get_breaks(reshp_NUTS_EU$y, plot_code)
palette = get_palette(plot_code, breaks)

# Map for Europe
eur_bbox = c(18.5,14.5,59.5,54.0)
eur_bbox = c(24,14.5,74,54) #_NUTS
europe = get_tm_shape(shp_country_all, eur_bbox, 
                      reshp_NUTS_EU, breaks, palette, alpha=1, midpoint=NA, legend_title="", legend=TRUE, 
                      shp_country_EU, title="TITLE", titlesize=1.0, titlefontface="bold", legend_size=0.6, 
                      frame=TRUE)

# Map for other regions
other_regions = list("azores"   = list("bbox"=c(9,22,14,28.5),     "title"="Azores (Portugal)"  ),
                     "canaries" = list("bbox"=c(15,9,20.75,12),    "title"="Canaries (Spain)"   ),
                     "cyprus"   = list("bbox"=c(62.5,15.5,66,18),  "title"="Cyprus"             )#,
                     # "maderia"  = list("bbox"=c(17,14.5,19.75,16), "title"="Maderia (Portugal)" ) # NUTS only
                     )
n_otherreg = length(other_regions)
for (or in 1:n_otherreg){
  other_regions[[or]]$map = get_tm_shape(shp_country_all, other_regions[[or]]$bbox,
                                         reshp_NUTS_EU, breaks, palette, alpha=NA, midpoint=NULL, legend_title=NA, legend=FALSE,
                                         shp_country_EU, other_regions[[or]]$title, titlesize=0.4, titlefontface="plain", legend_size=0.6, 
                                         frame=TRUE)
  }
other_regions[[or]]$map

# Export of the Map
pdf( paste0( foldout, "map_", plot_code, "_postmeta.pdf" ), width = 5, height = 5 )
print( europe )
print( other_regions$azores$map,   vp = viewport( 0.13, 0.76, width = 0.20, height = 0.20 ) )
print( other_regions$canaries$map, vp = viewport( 0.13, 0.60, width = 0.20, height = 0.10 ) )
print( other_regions$cyprus$map,   vp = viewport( 0.13, 0.49, width = 0.15, height = 0.09 ) )
dev.off()


final_shape =  tm_shape( shp_country_all, bbox=eur_bbox*10e4 ) +
  tm_fill( "#E0E0E0" ) + 
  tm_crs("auto") +
  tm_shape( shape_regions ) +
  tm_fill( "y", breaks=breaks, palette=palette, alpha=alpha, midpoint=midpoint, title=legend_title, legend.show=legend ) +
  tm_shape( shape_countries_EU ) +
  tm_borders( lwd = .25 ) +
  tm_layout( main.title=title, main.title.position=c("centre","top"), 
             main.title.size=titlesize, main.title.fontface=titlefontface, 
             legend.position=c(0.025,0.075), legend.text.size=legend_size, frame=TRUE )
c(24,14.5,74,54)

tm_shape(World, 
    bbox = sf::st_bbox(c(xmin = 24, xmax = 14.5, ymin = 74, ymax = 54))) +
  tm_polygons()