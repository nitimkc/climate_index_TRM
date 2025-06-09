# ------------------------------------------------------------------------------
# -
# - Author: Niti Mishra
# - Wrapper functions for obtaining data specific to phases of climate index,
#   running statistical tests
# - adapted from 
# -   NM as of April 5 2025 
# ------------------------------------------------------------------------------


# if (!require("pacman")) install.packages("pacman")
# suppressMessages(  pacman::p_load(config) )

# source("??");


get_detrended <- function(data, fit_form, y, x, loess_span=NA){
  
  na_rows = is.na(data[[y]])
  N = length(na_rows) - sum(na_rows)
  
  if (!is.na(loess_span)) {
    model = loess( fit_form, data[!na_rows,], span=loess_span/N ) 
  } else {
    model = lm( fit_form, data[!na_rows,] )
  }
  
  detrended = unlist( data[!na_rows, y] - predict(model, data[!na_rows, ]) )
  if (sum(na_rows)>0){
    detrended = insert( detrended, ats=which(na_rows), values=NA)
  }
  binary_threshold = detrended > 0 
  
  return(list("detrended"=detrended,
              "binary_threshold"=binary_threshold))
}

# get dates for specified period
get_period_dates <- function(all_dates, period="winter"){
  
  months = format(all_dates, "%m")
  if (period=="winter") {
    dates_period = all_dates[which((months=="12") | (months=="01") | (months=="02"))]
  } else if (period=="summer") {
    dates_period = all_dates[which((months=="07") | (months=="08"))]
  } else { stop("Invalid period type provided.")}
  
  print(length(dates_period))
  return(dates_period)
}


# data = attr_pos
# temp_group = "Total Heat"
get_group_mean <- function(temp_group, data, n_simu=1000){ 
  if (nrow(data)>0) {
    group_mean = data %>%
      filter( data[[temp_group]]==TRUE ) %>%
      select( 1:(n_simu+1) ) %>%
      colMeans(na.rm=TRUE)
  } else {
    group_mean = rep(NA, n_simu+1)
  }
  return(group_mean)
}


get_pval <- function(x, y, test_type="Student-t", paired=TRUE, alternative="two.sided"){
  if (! (any(is.na(x)) |  any(is.na(y))) ) {
    if (test_type=="Student-t"){
      pval = t.test(x, y, var.equal=FALSE, paired=paired, alternative=alternative)$p.value
    } else if (test_type=="Wilcox") {
      pval = wilcox.test(x, y, paired=paired, alternative=alternative)$p.value
    }
  } else {pval = NA}
  return(pval)
}
# ttest_stat = apply( test, 3, function(x) get_pval(x[1,], x[2,] ) )