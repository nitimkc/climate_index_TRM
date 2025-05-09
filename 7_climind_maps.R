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

# load required outputs from previous runs
fname_note = tail( strsplit(sub("\\.[^.]*$", "", CONFIG$NAO_THRESHOLD), "_")[[1]], 1)
test_type = "ttest"

# AF_NAO = readRDS( paste0(foldout, "AF_NAO_", fname_note, ".rds") )
# AF_NAO = readRDS( paste0(foldout, "AF_NAO_", fname_note, "_periods_parallelized.rds") )
# pval = readRDS( paste0(foldout, test_type, "_pval_AF_NAO_", fname_note, "_periods.rds") )

info_region  = readRDS( paste0(foldout, "info_region.rds") )
AF_NAO = readRDS( paste0(foldout, "AF_NAO_", fname_note, "_confint_phases_seasons.rds") )
pval = readRDS( paste0(foldout, test_type, "_pval_AF_NAO_", fname_note, "_seasons.rds"))

sig = pval <= 0.01
print(colSums(sig, na.rm=TRUE))
colnames(sig) = paste0("sig ", colnames(sig))

# others
temp_groups = strsplit(CONFIG$TEMP_RANGE, ",")[[1]]
n_groups  = length(temp_groups)
n_regions = length(info_region$code)
n_countries = length( unique(info_region$country_code) )
n_EU_regions = length( unique(info_region$EU_regions) )

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
# MAP ATTRIBUTABLE FRACTION
# ------------------------------------------------------------------------------

# plot vars
val = 1 # attr
geo_idx = n_countries+n_EU_regions
season = "summer"      # season
ssn_idx = ifelse(season=="winter", 1, ifelse(season=="summer", 2, 3) )


# Total temperature group only
# ----------------------------
start_time = Sys.time() # ~1 mins 10.59984 secs
grp = 1                # Total temp group
attr = AF_NAO[-seq_len(geo_idx),val,grp,ssn_idx,]
attr_sig = sig[,grp,ssn_idx]

plot_df = plot_df = as.data.frame( cbind( attr, "diff"=attr[,2]-attr[,1] ) )  # neg means positive phase has higher AF # how to handle NA
plot_df["sig"] = ifelse(attr_sig==T, "*", "")

for (phs in 0:1) { 
  
  if (phs==0){            
    fill_var = "diff"       # difference of phase or specific phase
    phase_name = paste0("_", fill_var)
    breaks  = cbreaks(range(plot_df[, fill_var], na.rm=TRUE), breaks_pretty(10))$breaks
    palette = "brewer.rd_bu" #-RdBu" # get_palette(breaks, plot_code="AFTH") # brewer.blues
    sig_var = NA #"sig"
    # eur_title = "Attributable Fraction between Positive and Negative phases of NAO"
    
    reshp_NUTS = base::merge(shp_NUTS, plot_df, by.x="NUTS_ID", by.y=0)      # add variables to plot
    eur = get_map(basemap, reshp_NUTS, fill_var, breaks, palette, midpoint=NA, alpha=0.7,
                  legend=TRUE, legend_title="", legend_size=0.4, frame=TRUE, sig=sig_var)
    insets = list()
    for (or in 1:length(other_regions)){
      print(other_regions[[or]]$title)
      reshp_NUTS_or = reshp_NUTS[(reshp_NUTS$NUTS_ID %in% other_regions[[or]]$islands), ]
      insets[[or]] = get_map(other_regions[[or]]$basemap, reshp_NUTS_or, fill_var, breaks, 
                             palette, midpoint=NA, alpha=0.7, legend=FALSE, frame=TRUE, 
                             title=other_regions[[or]]$title, sig=sig_var )
    }
    vp_list = list(viewport(x=0.12, y=0.80, width=0.15, height=0.15),
                   viewport(x=0.14, y=0.65, width=0.20, height=0.20),
                   viewport(x=0.10, y=0.45, width=0.20, height=0.20),
                   viewport(x=0.10, y=0.28, width=0.10, height=0.10))
    
    plotfname = paste0(foldout_plot, "AF_NAO_", fname_note, "_Total_", season, phase_name, ".png")
    tmap_save(eur, insets_tm=insets, insets_vp=vp_list, filename=plotfname, dpi=600)
    
  } else {
    fill_var = c("pos", "neg")    # each phase
    phase_name = "_phases"
    breaks  = cbreaks(range(plot_df[, fill_var], na.rm=TRUE), breaks_pretty(10))$breaks
    palette = "brewer.reds" # get_palette(breaks, plot_code="AFTH") # brewer.blues
    sig_var = NA
    # eur_title = "Attributable Fraction between Positive and Negative phases of NAO"
    
    reshp_NUTS = base::merge(shp_NUTS, plot_df, by.x="NUTS_ID", by.y=0)      # add variables to plot
    eur = get_map(basemap, reshp_NUTS, fill_var, breaks, palette, midpoint=NA, alpha=0.7,
                  legend=TRUE, legend_title="", legend_size=0.4, frame=TRUE, sig=sig_var)
    plotfname = paste0(foldout_plot, "AF_NAO_", fname_note, "_Total_", season, phase_name, ".png")
    tmap_save(eur, filename=plotfname, dpi=600)
  }
}
end_time = Sys.time()
print(end_time-start_time)


# All temperature groups
# ----------------------

fill_var = c("Total Heat", "Moderate Heat", "Extreme Heat", "Total Cold", "Moderate Cold", "Extreme Cold", "Total")
grps  = 1: length(temp_groups)       # all temp groups
attr = AF_NAO[-seq_len(geo_idx),val,grps,ssn_idx,]

for (phs in 0:2) { 
  start_time = Sys.time()   # ~2 mins
  
  if (phs==0){            # difference of phase or specific phase
    phase_name = "_diff" 
    plot_df = as.data.frame( attr[,,2] - attr[,,1] ) # neg means positive phase has higher AF # how to handle NA
    breaks  = cbreaks(range(plot_df, na.rm=TRUE), breaks_pretty(10))$breaks
    palette = "-RdBu"
    plot_df[, colnames(sig)[grps]] = ifelse(sig[,grps,season]==T, "*", "")
    plot_df[, colnames(sig)[grps]][is.na(plot_df[, colnames(sig)[grps]])] = ""
    sig_var = paste0("sig ", fill_var)
    
  } else {
    phase_name = ifelse(phs==1, "_pos", "_neg")
    plot_df = as.data.frame( attr[,,phs] )
    breaks  = cbreaks(range(plot_df, na.rm=TRUE), breaks_pretty(10))$breaks
    palette = "brewer.reds" # get_palette(breaks, plot_code="AFTH") # brewer.blues
    sig_var = NA
  }

  print(colMeans(plot_df[,grps], na.rm=T))
  reshp_NUTS = base::merge(shp_NUTS, plot_df, by.x="NUTS_ID", by.y=0)      # add variables to plot
  # eur_title = "Attributable Fraction between Positive and Negative phases of NAO"

  eur = get_map(basemap, reshp_NUTS, fill_var, breaks, palette=palette, midpoint=NA, alpha=0.7,
                legend=TRUE, legend_title="", legend_size=0.4, frame=TRUE, sig=sig_var)
  plotfname = paste0(foldout_plot, "AF_NAO_", fname_note, "_Groups_", season, phase_name, ".png")
  print(plotfname)
  tmap_save(eur, filename=plotfname, height = 8.27, width = 11.69, dpi=600)

  end_time = Sys.time()
  print(end_time-start_time)
}


# ------------------------------------------------------------------------------
# POPULATION WEIGHTED AF
# ------------------------------------------------------------------------------

# popn = read.table( paste0(CONFIG$ROOT, "indata/", "eurostat_table_popu_allsex_allage.csv"), 
#                    header = TRUE, sep = "," )
# popn = popn[popn$year==2024, c(1,3)]
# 
# grps = 1: length(temp_groups)
# n_simu = 1000
# attr = AF_NAO[,grps,2,1:n_simu] - AF_NAO[,grps,1,1:n_simu] # simu only
# 
# plot_df = as.data.frame(attr[,1,1])
# colnames(plot_df) = c("Total")
# plot_df = base::merge(popn, plot_df, by.x="location", by.y="row.names")
# plot_df["CNTR_CODE"] = substr(plot_df$location, 1, 2)
# 
# test = plot_df %>%
#   group_by(CNTR_CODE) %>%
#   summarise( Total = sum(popu * Total) / sum(popu) ) 
# 
# plot_df[, colnames(sig)[grps]] = ifelse(sig[,grps]==T, "*", "")
# plot_df[, colnames(sig)[grps]][is.na(plot_df[, colnames(sig)[grps]])] = ""
# plot_df["CNTR_CODE"] = substr(rownames(plot_df), 1, 2)
