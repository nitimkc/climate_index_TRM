# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra
# - Wrapper functions for generating NUTS EU maps
# - adapted from 
# -   Joan Ballester as of October 2024
# -   Tomás Janos as of February 25 2025
# -   NM as of April 2 2025 
# ------------------------------------------------------------------------------


if (!require("pacman")) install.packages("pacman")
suppressMessages( pacman::p_load(config, giscoR, lwgeom, sf, tmap, tmaptools) )

# source("??");


get_gisco_codes <- function(res, gisco, info){
  
  if( res=="region" | res=="country" ){ 
    gisco_revised = gisco[ gisco$NUTS_ID %in% info$code, ]
    
  } else if (grepl("NUT", res)){ # TO DO change to if "NUT" in res and remove NUTS list
    code_length = as.numeric( substr(res, nchar(res), nchar(res)) ) + 2 # last character + 2 
    gisco_revised = gisco[ gisco$NUTS_ID %in% info$code[which(nchar(info$code) == code_length)], ]
    
  } else { stop("Invalid case for resolution !!!"); }
  
  # gisco_revised = st_transform(gisco_revised)
  return( gisco_revised )
}
# gisco_plot = get_gisco_codes("region_NUT3", shp_NUTS_EU, info_NUTS) # 572 AT111
# gisco_plot = get_gisco_codes("region_NUT2", shp_NUTS_EU, info_NUTS) # 38  BG31
# gisco_plot = get_gisco_codes("region_NUT1", shp_NUTS_EU, info_NUTS) # 31  DE1
# gisco_plot = get_gisco_codes("country_NUT0", shp_NUTS_EU, info_NUTS)# 5   CY

get_breaks <- function(legends, constrained=NA, plot_code=NA, fraction=0.975 ){
  
  "plot_code: type of plot by code name
   fraction: Fraction of Values to be Represented in the Colorbar
  "
  if (is.na(constrained)){
    
    if ( grepl("RR", plot_code) ){ # Scale: Greater than One
      to   = ceiling( quantile(legends, fraction, na.rm=TRUE) * 10 ) / 10
      from = 1.000
      
    } else if ( grepl("AF", plot_code) ) { # AF Scale: Greater than Zero
      to   = ceiling( quantile(legends, fraction, na.rm=TRUE) )
      from = 0.000
    } else { stop("Invalid plot code!!!")}
    
  } else if (constrained==FALSE){ # unconstrained/ ?? NOT DIVERGENT 
    to   = ceiling( quantile(legends, 1 - (1-fraction)/2, na.rm=TRUE) )
    from =   floor( quantile(legends,     (1-fraction)/2, na.rm=TRUE) )
      
  } else if (constrained==TRUE){ # Scale: Symmetric Around Zero
    to   = ceiling( quantile(abs(legends), fraction, na.rm=TRUE) )
    from = -to
  }
  
  # interval step
  bar_range = to - from ; base10 = 10^( ceiling(log10(bar_range)) )
  if     ( base10 / 10 < bar_range & bar_range <= base10 / 5 ){ 
    by = base10 / 50 
  } else if( base10 /  5 < bar_range & bar_range <= base10 / 2 ){ 
    by = base10 / 20 
  } else if( base10 /  2 < bar_range & bar_range <= base10 / 1 ){ 
    by = base10 / 10 
  } else  { stop("Invalid by value for breaks !!!") }
  
  # Rounding the Colorbar Range Values According to the Colorbar Interval Value
  to = ceiling( to/by )*by; from = floor( from/by )*by
  return ( seq(from=from, to=to, by=by) )
}
# breaks = get_breaks(AF_NAO, plot_code="AFTH")


get_palette <-function(breaks, constrained=NA, plot_code=NA){
  " TO DO : CHANGE SHORT FORMS OF reg_plot to full info using paste0 and | and list reg_plots
  "
  N = length(breaks) - 1
  if (is.na(constrained)) {
    
    if( grepl("RR", plot_code) & as.numeric(gsub("\\D", "", plot_code))==75 ){ # White-to-Black Scale
      palette = colorRampPalette( c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252","#252525","#000000") )( N)
      
    } else if( grepl("RR", plot_code) & as.numeric(gsub("\\D", "", plot_code))<75 |  grepl("AFTC", plot_code)){ # White-to-Blue Scale
      palette = colorRampPalette( c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b") )( N )
      
    } else if( grepl("MM|AF|TA|TW|TS|TP|TI", plot_code) | grepl("RR", plot_code) & as.numeric(gsub("\\D", "", plot_code))>75 ){# White-to-Red Scale
      palette = colorRampPalette( c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d") )( N)
      
    } else { stop("Invalid plot code!!!")}
  
  } else if ( constrained==FALSE ){ # Blue-to-White-to-Red Scale
    palette = colorRampPalette( c("#08306b","#08519c","#2171b5","#4292c6","#6baed6","#9ecae1","#c6dbef","#deebf7","#f7fbff",
                                  "#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d") )( N)
  } else {"Invalid constrained value"}
  
  return(palette)
}
# palette = get_palette(breaks, plot_code="AFCH")


vec_to_bbox <- function(bbox, crs=3035){
  bbox = st_bbox( c(xmin=bbox[1],
                    ymin=bbox[2], 
                    xmax=bbox[3], 
                    ymax=bbox[4]), crs=crs )
  # bbox = st_as_sfc(bbox) # not required
  return( bbox )
}

# eur_bbox = vec_to_bbox( c(18.5,14.5,59.5,54.0)*10e4 )
# eur_bbox = vec_to_bbox( c(24,14.5,74,54)*10e4 ) 

get_basemap <- function(fill_shape, bbox, border_shape){ # add frame if needed
  
  basemap = tm_shape( fill_shape, bbox=bbox) +
    tm_fill( "#E0E0E0" ) +
    tm_shape( border_shape, is.main=FALSE ) +
    tm_borders( lwd = 1 )
  return(basemap)
}
# get_basemap(shp_country_all, or_bbox)


get_map <- function(basemap, fill_shape, fill_var, breaks, palette, midpoint=NULL, alpha=NA, 
                    legend=FALSE, legend_title=NA, legend_size=NA, frame=TRUE, title=NA, sig=NA) {
  
  final_shape = basemap  +
    tm_shape( fill_shape, is.main=FALSE) +
    tm_polygons( fill = fill_var, lwd = 0.5, lty = "blank", fill_alpha = alpha, fill.free = FALSE,
                 fill.scale = tm_scale_intervals(breaks=breaks, values=palette, midpoint=midpoint),
                 # fill.scale = tm_scale_intervals( values=palette, midpoint=NA),
                 fill.legend = tm_legend(show=legend, title=legend_title, title.size=legend_size), # legend.reverse=TRUE
                 ) 
    tm_layout( frame=frame )

  if ( !all(is.na(sig)) ){
    final_shape = final_shape +
      tm_text( text = sig, # currently uses the same col for all maps when using with facet
               size = .45,
               size.legend = tm_legend(show=FALSE),
               # size.scale = tm_scale_continuous(values.scale = 0.5)
               )
  }

  if (!is.na(title)){
    final_shape = final_shape +
      tm_title( title, size=0.4, fontface="bold")
  }
  
  return(final_shape)
}
