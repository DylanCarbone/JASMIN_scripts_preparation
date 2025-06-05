#' Wrapper function to fit occupancy models to multiple species.
#'#' Function to fit the occupancy models
#'
#'@param splist Vector of target species to model
#'@param obdata Data frame containing species occurrence records with the following columns: species, Date, gridref, Year, Week, Month, and optionally covnames (see below) and listL
#'@param occformula Formula for occupancy probability
#'@param detformula Formula for detection probability
#'@param covnames Vector of covariate names in obdata
#'@param seed Option random seed value to set
#'@param minyear Optional filter to first year of interest
#'@param maxyear Optional filter to last year of interest
#'@param trendyears Vector of start years for trend estimation. If \code{trendyears = NULL} then no trends will be calculated.
#'@param nstart Number of starting values to run. Default \code{nstart = 3}.
#'@param nrun Optional number of runs (since starting values are based on preceding year). Default \code{nrun = NULL} produces one run.
#'@param allsites Optional data frame of sites for which the occupancy index will be calculated for.
#'@param qval Quantile value to filter records to months where the species was observed. Default \code{qval = 0.025}.
#'@param prev_start Provide starting values e.g. based on outputs of a previous run.
#'@param outputdir Optional directory to save output files to.
#'@param printprogress print the progress of the run (only available for non-parallel option)
#'@param engine Choose the engine used by unmarked.
#'@return A list containing various outputs
#'@import data.table
#'@export

fit_occ_ms <- function(splist,
                      obdata,
                      occformula = "~North+I(North^2)+East+I(East^2)",
                      detformula = "~logLL+SEAS",
                      covnames = c("East","North"),
                      seed = NULL,
                      minyear = NULL,
                      maxyear = NULL,
                      trendyears = NULL,
                      nstart = 3,
                      nrun = 1,
                      allsites = NULL,
                      qval = 0.025,
                      prev_start = NULL,
                      outputdir = NULL,
                      printprogress = FALSE,
                      engine = "C"){

      if(!is.null(seed)) set.seed(seed)
      outputp <- list()
      for(ispp in splist){
        cat("Starting ", ispp," at ", base::date(),"\n")
          outputp[[ispp]] <- nrun_wrapper(ispp,
                                          nrun,
                                          obdata,
                                          occformula = occformula,
                                          detformula = detformula,
                                          covnames = covnames,
                                          minyear = minyear,
                                          maxyear = maxyear,
                                          trendyears = trendyears,
                                          allsites = allsites,
                                          qval = qval,
                                          nstart = nstart,
                                          printprogress = printprogress,
                                          prev_start = prev_start,
                                          engine = engine)

          cat("Finishing ", ispp," at ", base::date(),"\n")
        }
  return(outputp)
}

# Wrapper function to perform multiple model runs
nrun_wrapper <- function(spp,
                         nrun,
                         ...){
  if(nrun > 1){
    output <- list()
    for(irun in 1:nrun){
      cat("Starting ", spp, " run ", irun, " at ", base::date(), "\n")
      output[[irun]] <- fit_occ(spp,
                                irun = irun,
                                ...)
      }
  } else {
    output <- fit_occ(spp,
                      ...)
  }
  return(output)
}

