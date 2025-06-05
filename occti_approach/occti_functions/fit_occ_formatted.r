#' Function to fit the occupancy models
#'
#'@param spp Target species to estimate occupancy for.
#'@param obdata Data frame containing species occurrence records with the following columns: species, Date, gridref, Year, Week, Month, and optionally covnames (see below) and listL
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
#'@param printprogress Print the progress of the run (only available for non-parallel option)
#'@param engine Choose the engine used by unmarked.
#'@param prev_output Previous output to append.
#'@return A list containing various outputs
#'@import data.table
#'@import unmarked
#'@export

fit_occ_formatted <- function(spp,
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
                    printprogress = FALSE,
                    engine = engine,
                    prev_output = NULL,
                    irun = 0){

    # --- Initial setup ---
    # These variables are defined to avoid CRAN check notes for non-standard evaluation
    Year <- species <- Week <- N <- NULL

    # --- Set minimum and maximum year from data if not provided ---
    if (is.null(minyear)) {
    minyear <- min(obdata$Year)
    }

    if (is.null(maxyear)) {
    maxyear <- max(obdata$Year)
    }

    # --- Add list length (proxy for detection effort) if not already present ---
    if (!("listL" %in% colnames(obdata))) {
    obdata <- add_listL(obdata)
    }

    # --- Create list of all sites (with spatial covariates) if not provided ---
    if (is.null(allsites)) {
    allsites <- unique(obdata[, c("gridref", covnames), with = FALSE])
    }

    # Record start time for performance tracking
    st1 <- Sys.time()

    # --- Filter data to selected year range ---
    obdata <- obdata %>%
    filter(Year %in% minyear:maxyear)

    # --- Calculate seasonal detection probability by week for the species ---
    # Filter only rows for the target species
    obdata_sp <- obdata %>%
    filter(species == spp)

    # Count number of observations per week
    pweek <- obdata_sp %>%
    group_by(Week) %>%
    summarise(N = n(), .groups = "drop")

    # Convert weekly counts into proportions
    pweek <- pweek %>%
    mutate(N = N / sum(N))

    # If not all weeks are represented, add missing weeks with zero detection
    missing_weeks <- setdiff(unique(obdata$Week), pweek$Week)
    if (length(missing_weeks) > 0) {
    pweek <- bind_rows(
        pweek,
        data.frame(Week = missing_weeks, N = 0)
    )
    }

    # --- Combine weeks 52 and 53 ---
    # Assign value from week 52 to both week 52 and 53 for consistency
    pweek <- pweek %>%
    mutate(N = if_else(Week %in% c(52, 53),
                        pweek$N[pweek$Week == 52],
                        N))

    # --- Loop over years to fit model per year ---
    years <- obdata %>%
    filter(species == spp) %>%
    pull(Year) %>%
    unique() %>%
    sort()

    # Initialize output containers
    months <- coefs <- z <- aics <- NULL

    for (kyear in rev(years)) {
    if (printprogress) {
        cat(spp, "for", kyear, "at", base::date(), "\n")
    }

    # Subset site list for current year if "Year" column exists
    if ("Year" %in% colnames(allsites)) {
        allsitesk <- allsites %>%
        filter(Year == kyear)
    } else {
        allsitesk <- allsites
    }

    # Subset observation data for the current year
    obdatak <- obdata %>%
        filter(Year == kyear)

    # Skip loop iteration if there are no records for the species this year
    if (nrow(filter(obdatak, species == spp)) == 0) {
        message("no records for", kyear, "in species,", spp)
        next()
    }

    # --- Estimate active season based on observation months ---
    spp_months <- obdatak %>%
        filter(species == spp) %>%
        pull(Month)

    month1 <- if (is.null(qval)) {
        min(spp_months)
    } else {
        floor(quantile(spp_months, qval))
    }

    month2 <- if (is.null(qval)) {
        max(spp_months)
    } else {
        ceiling(quantile(spp_months, 1 - qval))
    }

    # Record the estimated active months for the year
    months <- bind_rows(months, data.frame(Year = kyear, min = month1, max = month2))

    # --- Filter records to active months and merge seasonal detection weights ---
    obdatak <- obdatak %>%
        filter(Month %in% month1:month2) %>%
        left_join(pweek, by = "Week")

    # --- Set up binary occurrence column: 1 if species present, 0 otherwise ---
    obdatak <- obdatak %>%
    mutate(Occ = if_else(species == spp, 1, 0))

    # --- Collapse to site-date visits, preserving max occurrence and unique effort ---
    obdatak1 <- obdatak %>%
    group_by(date, gridref, across(all_of(covnames)), Week, N) %>%
    summarise(
        Occ = max(Occ),
        listL = unique(listL),
        .groups = "drop"
    ) %>%
    arrange(gridref, desc(Occ))

    # --- Create site-by-visit matrices for occurrence, effort, and detection seasonality ---
    # These operations cannot be fully replaced by `dplyr`, so we retain base/paste + `plyr`
    obdatak1t <- do.call(
    plyr::rbind.fill.matrix,
    plyr::dlply(obdatak1, "gridref", function(a) { matrix(a$Occ, nrow = 1) })
    )

    obdatak1tL <- do.call(
    plyr::rbind.fill.matrix,
    plyr::dlply(obdatak1, "gridref", function(a) { matrix(a$listL, nrow = 1) })
    )

    obdatak1tPw <- do.call(
    plyr::rbind.fill.matrix,
    plyr::dlply(obdatak1, "gridref", function(a) { matrix(a$N, nrow = 1) })
    )

    # --- Extract spatial covariates per site ---
    obdatak1tEN <- obdatak1 %>%
    select(gridref, all_of(covnames)) %>%
    distinct()

    obdatak1tEN$gridref <- as.factor(obdatak1tEN$gridref)

    # --- Limit maximum number of visits per site to 50 ---
    if (ncol(obdatak1t) > 50) {
    obdatak1t   <- obdatak1t[, 1:50]
    obdatak1tL  <- obdatak1tL[, 1:50]
    obdatak1tPw <- obdatak1tPw[, 1:50]
    }

    # --- Create unmarkedFrame object for occupancy modeling ---
    dataf <- unmarkedFrameOccu(
    obdatak1t,
    obsCovs = list(
        logLL = log(obdatak1tL),
        SEAS  = obdatak1tPw
    ),
    siteCovs = obdatak1tEN
    )

    # --- Prepare starting values for model fitting ---
    occfit <- starts <- list()

    # Estimate number of parameters from both model formulas
    nparam <- length(attr(terms(formula(occformula)), "term.labels")) +
            length(attr(terms(formula(detformula)), "term.labels")) + 2

    # Set initial values
    if (is.null(prev_start)) {
    starts[[1]] <- rep(0, nparam)
    } else {
    starts[[1]] <- prev_start
    }

    # --- First model fit attempt ---
    occfit[["f1"]] <- try(
    occu(
        formula(paste(detformula, occformula, sep = "")),
        starts = starts[[1]],
        data = dataf,
        control = list(maxit = 1000),
        engine = engine
    ),
    silent = TRUE
    )

    # --- Additional model fits with random starts, if requested ---
    if (nstart > 1) {
    for (istart in 2:nstart) {
        if (is.null(prev_start)) {
        starts[[istart]] <- runif(nparam, -2, 0.2)
        } else {
        if (istart == 2) {
            starts[[istart]] <- rep(0, nparam)
        } else {
            starts[[istart]] <- prev_start + runif(nparam, -1, 1) * prev_start * 0.2
        }
        }

        occfit[[paste0("f", istart)]] <- try(
        occu(
            formula(paste(detformula, occformula, sep = "")),
            starts = starts[[istart]],
            data = dataf,
            control = list(maxit = 1000),
            engine = engine
        ),
        silent = TRUE
        )
    }
    }

    # --- Evaluate model fits using AIC, and select the best model ---
    aicsk <- rep(NA, length(occfit))

    for (i in seq_along(occfit)) {
    fit <- occfit[[i]]

    # Check if model fit or variance-covariance matrix failed
    if (inherits(fit, "try-error") ||
        inherits(try(unmarked::vcov(fit), silent = TRUE), "try-error") ||
        min(diag(unmarked::vcov(fit))) < 0 ||
        min(eigen(unmarked::vcov(fit, type = "state"))$values) < 0) {
        aicsk[i] <- NA
    } else {
        aicsk[i] <- fit@AIC
    }
    }

    # Skip year if all fits failed
    if (all(is.na(aicsk))) {
    next()
    }

    # Identify and retain best-fitting model
    best <- which(aicsk == min(aicsk, na.rm = TRUE))[1]
    occfit <- occfit[[best]]
    beststarts <- starts[[best]]

    # Save AIC values for all fits for later inspection
    aicsk <- data.frame(Year = kyear, start = 1:nstart, AIC = aicsk)
    aics <- bind_rows(aics, aicsk)

    # --- Save output if the model fit is valid ---
    if (class(occfit)[1] != "try-error" &&
        min(diag(unmarked::vcov(occfit))) > 0) {
    
    # Update starting values with the coefficients from the fitted model
    prev_start <- unmarked::coef(occfit)

    # Ensure state variance-covariance matrix is positive definite
    if (min(eigen(unmarked::vcov(occfit, type = "state"))$values) >= 0) {

        # --- Predict occupancy and estimate uncertainty ---

        # Extract the formula for state (occupancy) predictions
        stateformula <- as.formula(paste("~", occfit@formula[3], sep = ""))

        # Create model matrix for prediction using site covariates
        X <- model.matrix(stateformula, model.frame(stateformula, allsitesk))

        # Predict occupancy on logit scale
        y <- as.vector(X %*% unmarked::coef(occfit, "state"))
        test <- plogis(y)  # Convert to probability scale

        # Compute gradient for delta method
        dgdy <- plogis(y) / (1 + exp(y))
        dBeta <- X * dgdy
        mean_dBeta <- colMeans(dBeta, na.rm = TRUE)

        # Estimate variance and standard deviation of mean occupancy
        psi_var <- mean_dBeta %*% unmarked::vcov(occfit, type = "state") %*% mean_dBeta
        psi_sd <- sqrt(psi_var)

        # Mean occupancy across all sites
        I <- mean(test, na.rm = TRUE)

        # Delta method for logit-scale confidence interval
        mean_dBeta2 <- (1 / I + 1 / (1 - I)) * mean_dBeta
        psi_var_logit <- mean_dBeta2 %*% unmarked::vcov(occfit, type = "state") %*% mean_dBeta2

        # Define logit function
        logit <- function(x) { log(x / (1 - x)) }

        # --- Create summary data for year ---
        z1 <- data.frame(
        Year = kyear,
        psi = plogis(unmarked::coef(occfit)[1]),  # Intercept-based occupancy estimate
        psiA = I,  # Average predicted occupancy
        psiA_L = plogis(logit(I) - 1.96 * sqrt(psi_var_logit)),  # Lower CI (logit-transformed)
        psiA_U = plogis(logit(I) + 1.96 * sqrt(psi_var_logit)),  # Upper CI (logit-transformed)
        psiA_Lunbounded = I - 1.96 * psi_sd,  # Lower CI (linear scale, may be <0 or >1)
        psiA_Uunbounded = I + 1.96 * psi_sd,  # Upper CI (linear scale)
        psiA_SD = psi_sd,                    # Standard deviation of occupancy (linear scale)
        psiA_SDb = sqrt(psi_var_logit),      # Standard deviation (logit scale)
        AIC = occfit@AIC,
        nRecords = nrow(subset(obdatak, species == spp)),
        nSquares = length(unique(subset(obdatak, species == spp)$gridref)),
        month_min = month1,
        month_max = month2,
        nstart = nstart,
        beststart = best,
        nstartNA = sum(is.na(aicsk$AIC)),
        naic = uniqueN(round(aicsk[!is.na(aicsk$AIC), ]$AIC, 1)),
        irun = irun
        )

        # Append summary for current year
        z <- rbind(z, z1)

        # --- Save coefficient estimates with standard errors ---
        coefs <- rbind(
        coefs,
        data.frame(
            Year = kyear,
            Coef = names(unmarked::coef(occfit)),
            Est = unmarked::coef(occfit),
            SE = SE(occfit),
            starts = beststarts
        )
        )
    }
    }

  }

    # --- Record end time ---
    et1 <- Sys.time()

    # --- Combine outputs into a results list ---
    if (is.null(prev_output)) {
    
    # No previous output: start a fresh results list
    results <- list(
        species = spp,
        OccModel = occformula,
        DetModel = detformula,
        Index = z,
        Coefs = coefs,
        pweek = pweek,
        Totaltime = et1 - st1,
        minyear = minyear,
        maxyear = maxyear,
        months = months,
        qval = qval,
        nstart = nstart,
        aics = aics
    )
    
    # Add run index if provided
    if (!is.null(irun)) {
        results$irun <- irun
    }

    } else {
    
    # Previous output exists: merge new results with existing ones
    z <- rbind(
        prev_output$Index[!prev_output$Index$Year %in% years, ],
        z
    )
    
    results <- list(
        species = spp,
        OccModel = occformula,
        DetModel = detformula,
        Index = z,
        Coefs = rbind(
        prev_output$Coefs[!prev_output$Coefs$Year %in% years, ],
        coefs
        ),
        pweek = pweek,
        Totaltime = et1 - st1,
        minyear = minyear,
        maxyear = maxyear,
        months = months,
        qval = qval,
        irun = NULL  # Reset irun if merging multiple runs
    )
    }

  # HERE

# --- Estimate temporal trends if requested ---
if (!is.null(trendyears) && !is.null(z) && nrow(z) > 2) {
  # Fit trend models for each year in trendyears (after minyear)
  trends <- do.call(
    rbind,
    lapply(trendyears[trendyears >= minyear], fit_trend, z = z, endyear = maxyear)
  )
} else {
  trends <- NULL
}

    # Add trend results to the final output
    results$Trends <- trends
    results$Trendyears <- trendyears

    # --- Return final results list ---
    return(results)
}
