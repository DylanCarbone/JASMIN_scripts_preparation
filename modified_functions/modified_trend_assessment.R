
modified_trend_assessment <- function (dat, summary, method = "lambda", start_year, end_year, 
                                       species_stat = "mean") 
  # LB added summary as an extra argument: BRCindicators:::indicator assessment seems to need it
{
  if (!method %in% c("lambda", "bma")) 
    stop("Method must be one of 'lambda' or 'bma'")
  if (method == "lambda") {
    #start_year <- min(dat$summary$year) LB removed: these are arguments to the function: overriding them means short term trends aren't estimated
    #end_year <- max(dat$summary$year)
    sp_assess <- BRCindicators:::species_assessment(dat = dat$LogLambda, 
                                    method = method, start_year = start_year, end_year = end_year, # LB replace start_year = start_year +1 with start_year =start_year
                                    species_stat = species_stat, plot = FALSE)
    ind_assessment <- BRCindicators:::indicator_assessment(summary_table = summary, 
                                           start_year = start_year, end_year = end_year)
    return(list(species_assessment = sp_assess, indicator_asssessment = ind_assessment)) 
  }
  else {
    #start_year <- min(dat$year) # LB removed: these are arguments to the function: overriding them means short term trends aren't estimated
    #end_year <- max(dat$year)
    names(dat) <-gsub("Year", "year", names(dat)) # the column name has to be year (lowercase) inside species_assessment
    
    sp_assess <- BRCindicators:::species_assessment(dat = dat, method = method, 
                                    start_year = start_year, end_year = end_year, species_stat = species_stat, 
                                    plot = FALSE)
    ind_assessment <- BRCindicators:::indicator_assessment(summary_table = summary, 
                                           start_year = start_year, end_year = end_year)
    return(list(species_assessment = sp_assess, indicator_asssessment = ind_assessment))
  }
}
