#' Function to estimate occupancy trends
#'@param startyear First year to calculate the trend from
#'@param z Data frame output from \code{fit_occ}, containing estimated occupancy index and standard errors.
#'@param endyear Optional. Year to calculate the trend up to. Default will select the last year of the index.
#'@param predictions Option to include predicted values in the output e.g. for plotting.
#'@export

# Estimate trends
fit_trend <- function(startyear, z, endyear=NULL, predictions = FALSE){

  if(quantile(z$psiA, .9, na.rm=TRUE) < .9){
    z <- subset(z, psiA < .99)
  }
  if(is.null(endyear)) endyear <- max(z$Year)
  z <- subset(z, Year %in% startyear:endyear)
  trd <- NULL
  fit3 <- NA
  if(nrow(z)>=3 && (min(z$Year)-startyear) <= 5){
    fit1 <- try(glm(psiA~Year, data=z, weights=1/psiA_SD, family=binomial(link="logit")), silent=TRUE)
     if(class(fit1)[1] != "try-error"){
      fit1.pred <- predict(fit1,
                           newdata=data.frame(Year=startyear:endyear),
                           weights=1/z$psiA_SD,
                           se.fit=TRUE)
      upp <- expit(fit1.pred$fit + 1.96 * fit1.pred$se.fit)
      pred <- expit(fit1.pred$fit)
      low <- expit(fit1.pred$fit - 1.96 * fit1.pred$se.fit)
      trd_ci1 <- round(pcfunc2(head(upp,1),tail(low,1)),3)
      trd_ci2 <- round(pcfunc2(head(low,1),tail(upp,1)),3)
      trd_ci <- sort(c(trd_ci1,trd_ci2))
      trd_se <- (trd_ci[2]-trd_ci[1])/1.96/2

      trd <- data.frame(Startyear=startyear,
                      Endyear=endyear,
                      nyears=nrow(z),
                      minyear=min(z$Year),
                      trd = round(pcfunc(pred),3),
                      trd_se = trd_se,
                      trd_ci1 = trd_ci1,
                      trd_ci2 = trd_ci2,
                      coef = coef(summary(fit1))[2,1],
                      SE = coef(summary(fit1))[2,2],
                      pval = coef(summary(fit1))[2,4])

    trd$sig <- if(!is.nan(trd$pval)){sigfunc(trd$pval)} else {NA}

    if(predictions){
      trd <- list(trd = trd,
                  predictions = data.frame(Year=startyear:endyear,
                                           pred = pred,
                                           low = low,
                                           upp = upp))
    }

  }}

  return(trd)
}

