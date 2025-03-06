#' Function to fit the occupancy models
#'
#'@param spp Target species to estimate occupancy for.
#'@param obdata Data frame containing species occurrence records with the following columns: Species, Date, Gridref, Year, Week, Month, and optionally covnames (see below) and listL
#'@param occformula Formula for occupancy probability
#'@param detformula Formula for detection probability
#'@param covnames Vector of covariate names in obdata
#'@param minyear First year of interest
#'@param maxyear Last year of interest
#'@param trendyears Vector of start years for trend estimation. If \code{trendyears = NULL} then no trends will be calculated.
#'@param nstart Number of starting values to run. Default \code{nstart = 3}.
#'@param allsites Optional data frame of sites for which the occupancy index will be calculated for.
#'@param qval Quantile value to filter records to months where the species was observed. Default \code{qval = 0.025}.
#'@param prev_start Provide starting values e.g. based on outputs of a previous run.
#'@param outputdir Optional directory to save output files to.
#'@param printprogress Print the progress of the run (only available for non-parallel option)
#'@param engine Choose the engine used by unmarked.
#'@param prev_output Previous output to append.
#'@return A list containing various outputs
#'@import data.table
#'@import unmarked
#'@export

# Fit occupancy models
fit_occ <- function(spp,
                    obdata,
                    occformula = "~North+I(North^2)+East+I(East^2)",
                    detformula = "~logLL+SEAS",
                    covnames = c("East","North"),
                    minyear = NULL,
                    maxyear = NULL,
                    trendyears = NULL,
                    nstart = 1,
                    allsites = NULL,
                    qval = NULL,
                    prev_start = NULL,
                    outputdir = NULL,
                    printprogress = FALSE,
                    engine = engine,
                    prev_output = NULL,
                    irun = 0){

  # Satisfy not finding global variable
  Year <- Species <- Week <- N <- NULL

  #obdata <- as.data.table(obdata)

  if(is.null(minyear)) minyear <- min(obdata$Year)
  if(is.null(maxyear)) maxyear <- max(obdata$Year)

  # Add a column for list length is needed
  if(!("listL" %in% colnames(obdata))) obdata <- add_listL(obdata)

  if(is.null(allsites)) # allsites <- unique(select(obdata, c("Gridref", covnames)))
    allsites <- unique(obdata[, c("Gridref", covnames), with=FALSE])

  st1 <- Sys.time()
  # Data prep
  #==================================================================
  #obdata <- filter(obdata, Year %in%  minyear:maxyear)
  obdata <- obdata[obdata$Year %in%  minyear:maxyear,]

  # Calculate seasonal variation across years
  obdata_sp <- filter(obdata, Species == spp)
  pweek <- obdata_sp %>% group_by(Week) %>% summarise(N=n())
  #obdata_sp <- obdata[Species == spp]
  #pweek <- obdata_sp[, .N, by = Week]
  pweek$N <- pweek$N/sum(pweek$N)
  if(nrow(pweek) < length(unique(obdata$Week))){
    pweek <- rbind(pweek, data.frame(Week =
                                       unique(obdata$Week)[!(unique(obdata$Week)
                                                             %in% pweek$Week)], N=0))
  }
  # Make weeks 52 and 53 the same
  pweek[pweek$Week %in% c(52,53),]$N <- pweek[pweek$Week == 52,]$N
  # pweek[Week %in% c(52,53), N := pweek[pweek$Week == 52,]$N]

  # Loop over years
  #==================================================================
  years <- sort(unique(obdata[obdata$Species == spp,]$Year))
  months <- coefs <- z <- aics <- NULL
  for(kyear in  rev(years)){
    m <- 0
    if(printprogress)cat(spp,"for",kyear,"at",base::date(),"\n")

    if("Year" %in% colnames(allsites)){
      allsitesk <- filter(allsites, Year == kyear)
      #allsitesk <- allsites[Year == kyear]
    } else {
      allsitesk <- allsites}

    # NMRS data for selected year
    obdatak <- filter(obdata, Year == kyear)
    #obdatak <- obdata[Year == kyear]

    # Check there are records that year!
    if(nrow(subset(obdatak, Species==spp))==0) next()

    # Limit to species months
    month1 <- if(is.null(qval)){
                    min(obdatak[obdatak$Species == spp,]$Month)
                   } else {
                    floor(quantile(obdatak[obdatak$Species == spp,]$Month, qval))
                     }
    month2 <- if(is.null(qval)){
                    max(obdatak[obdatak$Species == spp,]$Month)
                   } else {
                    ceiling(quantile(obdatak[obdatak$Species == spp,]$Month, 1 - qval))
                     }
  months <- rbind(months, data.frame(Year=kyear,
                                       min=month1,
                                       max=month2))
    obdatak <- obdatak[obdatak$Month %in% month1:month2,]
    obdatak <- merge(obdatak, pweek, by = "Week")

    # Data prep
    #==================================================================
    # Extract data for species of interest
    obdatak$Occ <- 0
    obdatak[obdatak$Species == spp,]$Occ <- 1

    # Add list length info
    obdatak1 <- obdatak[,list(Occ=max(Occ), listL=unique(listL)), by=c("Date", "Gridref", covnames, "Week", "N")]

    # Needs ordering appropriately by site then date I think
    obdatak1 <- obdatak1[order(obdatak1$Gridref, -obdatak1$Occ),]

    # Reduce to minimal number of visits
    obdatak1t <- do.call(plyr::rbind.fill.matrix, plyr::dlply(obdatak1,
                                                              "Gridref", function(a){matrix(a$Occ, nrow=1)}))
    obdatak1tL <- do.call(plyr::rbind.fill.matrix, plyr::dlply(obdatak1,
                                                               "Gridref", function(a){matrix(a$listL, nrow=1)}))
    obdatak1tPw <- do.call(plyr::rbind.fill.matrix, plyr::dlply(obdatak1,
                                                                "Gridref", function(a){matrix(a$N, nrow=1)}))
    obdatak1tEN <- unique(obdatak1[,c("Gridref", covnames), with=FALSE])

    obdatak1tEN$Gridref <- as.factor(obdatak1tEN$Gridref)

    # Reduce max no. visits to 50
    if(ncol(obdatak1t) > 50){
      obdatak1t <- obdatak1t[,1:50]
      obdatak1tL <- obdatak1tL[,1:50]
      obdatak1tPw <- obdatak1tPw[,1:50]
    }

    # Data set up
    dataf <- unmarkedFrameOccu(obdatak1t,
                               obsCovs = list(logLL=log(obdatak1tL), SEAS=obdatak1tPw),
                               siteCovs = obdatak1tEN)

    # Fit occupancy model
    # if(multistart){
    occfit <- starts <- list()
    nparam <- length(attr(terms(formula(occformula)),"term.labels"))+
    length(attr(terms(formula(detformula)),"term.labels"))+2
    if(is.null(prev_start)){
        # If starting values not provided then try zeros and two other random starts
             starts[[1]] <- rep(0, nparam)
          } else {
             starts[[1]] <- prev_start
          }
        occfit[["f1"]] <- try(occu(formula(paste(detformula, occformula, sep="")),
                                   starts = starts[[1]],
                                   dataf, control=list(maxit=1000), engine = engine), silent=TRUE)
      if(nstart > 1){
        for(istart in 2:nstart){
          if(is.null(prev_start)){
                starts[[istart]] <- runif(nparam, -2., .2)
              } else {
                # First test a random start
                if(istart == 2){
                  starts[[istart]] <- rep(0, nparam)
                } else {
                # Then slight alternatives to prev_start
                  starts[[istart]] <- prev_start + runif(nparam, -1, 1)*prev_start*.2
                }
              }
          occfit[[paste0("f", istart)]] <- try(occu(formula(paste(detformula, occformula, sep="")),
                                            starts = starts[[istart]],
                                           dataf, control=list(maxit=1000), engine = engine), silent=TRUE)
        }}

      aicsk <- rep(NA, length(occfit))
      for(i in 1:length(occfit)){
        if(class(occfit[[i]])[1] =="try-error" ||
           class(try(unmarked::vcov(occfit[[i]]), silent = TRUE))[1] == "try-error" ||
           min(diag(unmarked::vcov(occfit[[i]]))) < 0 ||
           min(eigen(unmarked::vcov(occfit[[i]], type="state"))$values) < 0){
          #occfit[[i]] <- NULL
          aicsk[i] <- NA
        } else {
          aicsk[i] <- occfit[[i]]@AIC
        }
      }
      # If null for all starts tried then skip this year
      if(all(is.na(aicsk))) next()
      # Save the best model in terms of aic
      best <- which(aicsk == min(aicsk, na.rm=TRUE)[1])
      occfit <- occfit[[best]]
      beststarts <- starts[[best]]

      aicsk <- data.frame(Year = kyear, start = 1:nstart, AIC = aicsk)
      aics <- rbind(aics, aicsk)


    # Error checking
    #==================================================================

    if(is.null(occfit)) next()
    if(class(occfit)[1] != "try-error" && min(diag(unmarked::vcov(occfit))) > 0){
      prev_start <- unmarked::coef(occfit)
      if(min(eigen(unmarked::vcov(occfit, type="state"))$values)>=0){
        # Estimate index and standard error
        #==================================================================
        stateformula <- as.formula(paste("~", occfit@formula[3], sep=""))
        X <- model.matrix(stateformula, model.frame(stateformula, allsitesk))
        y <-  c(X %*% unmarked::coef(occfit, "state"))
        test <- plogis(y)
        dgdy <- plogis(y)/(1+exp(y)) #  exp(y)/(1+exp(y))^2
        dBeta <- X*dgdy
        mean_dBeta <- colMeans(dBeta, na.rm=TRUE)
        psi_var <- mean_dBeta %*% unmarked::vcov(occfit, type="state") %*% mean_dBeta
        psi_sd <- sqrt(psi_var)


        I <- mean(test, na.rm=TRUE)
        mean_dBeta2 <- (1/I + 1/(1-I))*mean_dBeta
        psi_var_logit <- mean_dBeta2 %*% unmarked::vcov(occfit, type="state") %*% mean_dBeta2

        `logit` <- function(x){ log(x/(1-x)) }

        z1 <- data.frame(Year = kyear,
                         psi = plogis(unmarked::coef(occfit)[1]),
                         psiA = I,
                         psiA_L = plogis(logit(I) - 1.96*sqrt(psi_var_logit)),
                         psiA_U = plogis(logit(I) + 1.96*sqrt(psi_var_logit)),
                         psiA_Lunbounded = I - 1.96*psi_sd,
                         psiA_Uunbounded = I + 1.96*psi_sd,
                         psiA_SD = psi_sd,
                         psiA_SDb = sqrt(psi_var_logit),
                         AIC = occfit@AIC,
                         nRecords = nrow(subset(obdatak, Species == spp)),
                         nSquares = length(unique(subset(obdatak, Species == spp)$Gridref)),
                         month_min = month1,
                         month_max = month2,
                         nstart = nstart,
                         beststart=best,
                         nstartNA = sum(is.na(aicsk$AIC)),
                         naic = uniqueN(round(aicsk[!is.na(aicsk$AIC),]$AIC, 1)),
                         irun = irun)

        if(!is.null(irun))
            z1$irun <- irun

        z <- rbind(z, z1)

        coefs <- rbind(coefs, data.frame(Year = kyear,
                                         Coef = names(unmarked::coef(occfit)),
                                         Est = unmarked::coef(occfit),
                                         SE = SE(occfit),
                                         starts = beststarts))

        #setDT(coefs)
        #coefsm <- coefs[, .(Est = mean(Est)), by = Coef]
        #prev_start <- coefsm$Est
        #names(prev_start) <- coefsm$Coef
      }}
  }

  et1 <- Sys.time()

  if(is.null(prev_output)){
    results <- list(Species = spp,
                  OccModel = occformula,
                  DetModel = detformula,
                  Index = z,
                  Coefs = coefs,
                  pweek = pweek,
                  Totaltime = et1 - st1,
                  minyear = minyear,
                  maxyear = maxyear,
                  months = months,
                  qval=qval,
                  nstart = nstart,
                  aics = aics)
    if(!is.null(irun))
      results$irun <- irun
  } else {
    z <- rbind(prev_output$Index[!prev_output$Index$Year %in% years,], z)

    results <- list(Species = spp,
                    OccModel = occformula,
                    DetModel = detformula,
                    Index = z,
                    Coefs = rbind(prev_output$Coefs[!prev_output$Coefs$Year %in% years,], coefs),
                    pweek = pweek,
                    Totaltime = et1 - st1,
                    minyear = minyear,
                    maxyear = maxyear,
                    months = months,
                    qval=qval,
                    irun = NULL)

  }

  if(!is.null(trendyears) && !is.null(z) && nrow(z) > 2){
    trends <- do.call(rbind,
                      lapply(trendyears[trendyears >= minyear],
                             fit_trend, z = z, endyear = maxyear))
  } else {
    trends <- NULL}

  results$Trends <- trends
  results$Trendyears <- trendyears

  if(!is.null(outputdir))
    saveRDS(results,
            file=paste(outputdir,
                       spp,"_occupancy_output.rds",sep=""))

  return(results)
}





