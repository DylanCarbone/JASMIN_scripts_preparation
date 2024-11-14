modified_calcMSI <- function (dat, method, write, outPath, plotLabel, bmaInd = NULL, 
                              startYear = NULL, endYear = NULL) # LB added startYear and endYear arguments
  #otherwise indicator is not set to 100 i
{
  if (!method %in% c("lambda", "bma")) 
    stop("Method must be one of lambda or bma")
  year_cols <- grep("year_", 
                   colnames(dat))
  colnames(dat) <- gsub("year_", 
                        "", colnames(dat))
  minYr <- ifelse(is.null(startYear), 
                  yes = min(as.integer(colnames(dat)[year_cols])),
                  no = startYear)
  maxYr <- ifelse(is.null(endYear), 
                  yes = max(as.integer(colnames(dat)[year_cols])),
                  no = endYear)
    if (method == "lambda") {
    
    
    #arr <- wrappeR:::sampArray(dat = dat, startYear = minYr, endYear = maxYr) # sampArray doesn't remove years
    # the columns are named as years so can just index by those?
    arr <- reshape2::acast(melt(dat, id = c("iteration", "species")), 
                           species ~ variable ~ iteration, value.var = "value")
    #start <- (startYear - startYear) + 1
    #end <- (endYear - startYear) + 1
    arr <- arr[, as.character(minYr:maxYr), ] 
    ind <- BRCindicators::lambda_indicator(input = arr, index = 100, 
                                           threshold_sd = Inf, threshold_Rhat = Inf, threshold_yrs = 10, 
                                           upperQuantile = 0.975, lowerQuantile = 0.025)
    summary <- ind$summary
    lt <- modified_trend_assessment(ind, summary = summary, start_year = min(ind$summary$year), 
                                    end_year = max(ind$summary$year))
    st <- modified_trend_assessment(ind, summary = summary, start_year = (max(ind$summary$year) - 
                                                         5), end_year = max(ind$summary$year))
    final <- modified_trend_assessment(ind, summary = summary, start_year = max(ind$summary$year), 
                                       end_year = max(ind$summary$year))
  }
  else {
    fudgeOcc <- function(x, fudgeFac = 1e-04) {
      x[x == 0] <- fudgeFac
      x[x == 1] <- 1 - fudgeFac
      return(x)
    }
    dat <- dat[, c(year_cols, "species")]
    dat[, year_cols] <- apply(dat[, year_cols], 2, function(x) {
      log(fudgeOcc(x))
    })
    # insert year filter here
    means <- aggregate(. ~ species, data = dat, mean, na.action = na.pass)
    #table(reshape::melt(means)$species, useNA = "always")
    sds <- aggregate(. ~ species, data = dat, sd, na.action = na.pass)
    getSumStats <- function(stat) {
      out <- reshape2::melt(stat, id.vars = "species", 
                            variable.name = "year", value.name = "index")
      return(out)
    }
    means <- getSumStats(means)
    se <- getSumStats(sds)[, 3]
    inDat <- data.frame(means, se = se)
    inDat$se[inDat$se == 0] <- 1e-04
    inDat$year <- as.numeric(gsub("year_", "", 
                                  inDat$year))
    ind <- BRCindicators::bma(data = inDat, seFromData = TRUE, 
                              m.scale = "loge")
    if (bmaInd != "prime") {
      summary <- data.frame(indicator = ind$Index.M, 
                            lower = ind$lowerCI.M, 
                            upper = ind$upperCI.M, 
                            Species_Number = as.integer(table(inDat[complete.cases(inDat), "year"])), #LB replaced this with complete.cases(inDat) to reflect the number of species with an index estimate per year 
                            year = ind$Year) # also added year column as needed by BRCindicators:::indicator_assessment
    }
    else {
      summary <- data.frame(indicator = ind$Index.Mprime, 
                            lower = ind$lowerCI.Mprime, 
                            upper = ind$upperCI.Mprime, 
                            Species_Number = as.integer(table(inDat[complete.cases(inDat), "year"])), # should this be the number of species with data each year?
                            year = ind$Year)
    }
    lt <- modified_trend_assessment(dat = ind, summary = summary, method = method, # changed first argument from ind to summary
                                          start_year = min(ind$Year), end_year = max(ind$Year))
    st <- modified_trend_assessment(dat = ind, summary = summary, method = method, 
                                          start_year = (max(ind$Year) - 5), end_year = max(ind$Year))
    final <- modified_trend_assessment(dat = ind, summary = summary, method = method, 
                                             start_year = max(ind$Year), end_year = max(ind$Year))
  }
  out <- list(summary, ind, st, lt, final)
  names(out) <- c("Summary", "MetaData", "st", 
                  "lt", "final")
  return(out)
}

