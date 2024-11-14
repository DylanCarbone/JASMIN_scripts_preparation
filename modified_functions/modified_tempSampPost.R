# indata = paste0(roster$modPath, roster$group,
#                 "/occmod_outputs/", roster$ver, "/"); keep = keep;
#                 keep_iter = keep_iter; output_path = NULL; region = roster$region;
#                 sample_n = roster$nSamps; group_name = roster$group;
#                 combined_output = TRUE; write = FALSE; minObs = roster$minObs;
#                 scaleObs = roster$scaleObs; t0 = roster$t0; tn = roster$tn;
#                 parallel = parallel; filetype = filetype

modified_tempSampPost <- 
function (indata = "../data/model_runs/", keep, keep_iter, output_path = "../data/sampled_posterior_1000/", 
          region, sample_n = 999, tolerance = 0, group_name = "", combined_output = TRUE, 
          max_year_model = NULL, min_year_model = NULL, write = FALSE, 
          minObs = NULL, scaleObs = "global", t0, tn, parallel = TRUE, 
          n.cores = NULL, filetype = "rdata") 
{
  if (parallel & is.null(n.cores)) 
    n.cores <- parallel::detectCores() - 1
  REGION_IN_Q <- paste0("psi.fs.r_", region)
  # LB moved this section to within combSamps function, as min iter is sometimes different for some species (see Bryophytes) 
  # if (!is.null(keep_iter)) { 
    # findMinIteration <- function(list_of_file_names) {
    #   if (length(list_of_file_names) < 1)
    #     stop("Error: list_of_file_names is empty")
    #   if (!is.character(list_of_file_names))
    #     stop("Error: list_of_file_names must be a character")
    #   list_of_file_names <- gsub("_[[:digit:]]{1}$", "",
    #                              list_of_file_names)
    #   iterations <- regmatches(list_of_file_names, regexpr("[[:digit:]]+$",
    #                                                        list_of_file_names))
    #   return(min(as.numeric(iterations)))
    # }
    # min_iter <- findMinIteration(keep_iter)
  #}
  samp_post <- NULL
  load_rdata <- function(fileName) {
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  combineSamps <- function(species, minObs) {
    out_dat <- NULL
    out_meta <- NULL
    raw_occ <- NULL
    nRec_glob <- NA
    nRec_reg <- NA
    nRec <- NA
    rot <- NULL
    if (!is.null(keep_iter)) {
      # LB moved this here (see comment above)
      findMinIteration <- function(list_of_file_names) { 
        if (length(list_of_file_names) < 1)
          stop("Error: list_of_file_names is empty")
        if (!is.character(list_of_file_names))
          stop("Error: list_of_file_names must be a character")
        list_of_file_names <- gsub("_[[:digit:]]{1}$", "",
                                   list_of_file_names)
        iterations <- regmatches(list_of_file_names, regexpr("[[:digit:]]+$",
                                                             list_of_file_names))
        return(min(as.numeric(iterations)))
      }
      min_iter <- findMinIteration(grep(species, keep_iter, value = TRUE)) # modified to extract a unique min_iter for each species 
      
      try(out_dat <- load_rdata(paste0(indata, species, 
                                       "_20000_1.rdata")))
      try(out_meta <- load_rdata(paste0(indata, species, 
                                        "_", min_iter, "_1.rdata")))
    }
    else {
      if (filetype == "rds") {
        try(out_dat <- readRDS(paste0(indata, species, 
                                      ".rds")))
        out_meta <- out_dat
      }
      else if (filetype == "rdata") {
        try(out_dat <- load_rdata(paste0(indata, species, 
                                         ".rdata")))
        out_meta <- out_dat
      }
    }
    if (!is.null(out_dat$model) & !is.null(out_meta)) {
      nRec_glob <- out_meta$species_observations
      if (scaleObs != "global") {
        dat <- out_meta$model$data()
        if (region %in% out_meta$regions) {
          region_site <- dat[[paste0("r_", region)]][dat$Site]
        }
        else {
          region_aggs <- unlist(out_meta$region_aggs[region])
          region_site <- rowSums(sapply(region_aggs, 
                                        function(x) dat[[paste0("r_", x)]][dat$Site]))
        }
        dat <- data.frame(year = dat$Year, rec = dat$y, 
                          region_site = region_site)
        dat <- dat[dat$region_site == 1 & dat$year >= 
                     (t0 - (out_meta$min_year - 1)) & dat$year <= 
                     (tn - (out_meta$min_year - 1)), ]
        nRec_reg <- sum(dat$rec)
      }
    }
    if (scaleObs == "global") 
      nRec <- nRec_glob
    else nRec <- nRec_reg
    print(paste0("load: ", species, ", ", scaleObs, " records: ", 
                 nRec))
    if(is.null(out_meta$region_aggs) | length(out_meta$region_aggs) == 0){ # LB 
      REGION_IN_Q <- "psi.fs"      # LB: for UK-scale models with no region data 
    }
    # LB
    # LB replaced:
    # if (nRec >= minObs & region %in% c(out_meta$regions, 
    #                                    names(out_meta$region_aggs)) & !is.null(out_dat$model)) {
    # with:
    if (nRec >= minObs & !is.null(out_dat$model)) { #LB, from here onwards the same 
      if (!is.null(keep_iter)) {
        out_dat1 <- NULL
        out_dat2 <- NULL
        out_dat3 <- NULL
        try(out_dat1 <- load_rdata(paste0(indata, species, 
                                          "_20000_1.rdata")))
        raw_occ1 <- data.frame(out_dat1$BUGSoutput$sims.list[REGION_IN_Q])
        try(out_dat2 <- load_rdata(paste0(indata, species, 
                                          "_20000_2.rdata")))
        raw_occ2 <- data.frame(out_dat2$BUGSoutput$sims.list[REGION_IN_Q])
        try(out_dat3 <- load_rdata(paste0(indata, species, 
                                          "_20000_3.rdata")))
        raw_occ3 <- data.frame(out_dat3$BUGSoutput$sims.list[REGION_IN_Q])
        if (!is.null(out_dat1) & !is.null(out_dat2) & 
            !is.null(out_dat3)) 
          raw_occ <- rbind(raw_occ1, raw_occ2, raw_occ3)
      }
      else {
        raw_occ <- data.frame(out_dat$BUGSoutput$sims.list[REGION_IN_Q])
      }
      if (!is.null(raw_occ)) {
        if (!is.null(keep_iter)) {
          diff <- (out_dat$BUGSoutput$n.sims * 3) - sample_n
        }
        else {
          diff <- out_dat$BUGSoutput$n.sims - sample_n
        }
        if (diff > tolerance) {
          raw_occ <- raw_occ[sample(1:nrow(raw_occ), 
                                    sample_n), ]
        }
        else if (abs(diff) <= tolerance) {
          print(paste("no sampling required: n.sims =", 
                      out_dat$BUGSoutput$n.sims))
        }
        else stop("Error: Not enough iterations stored. Choose a smaller value of sample_n")
        colnames(raw_occ) <- paste("year_", out_meta$min_year:out_meta$max_year, 
                                   sep = "")
        raw_occ$iteration <- 1:sample_n
        raw_occ$species <- species
        if (combined_output != TRUE) {
          write.csv(raw_occ, file = paste(output_path, 
                                          gsub(".rdata", "", i), "_sample_", sample_n, 
                                          "_post_", REGION_IN_Q, ".csv", sep = ""), 
                    row.names = FALSE)
        }
        out1 <- raw_occ
        dat <- out_meta$model$data()
        dat <- data.frame(year = dat$Year, rec = dat$y)
        first <- min(dat$year[dat$rec == 1]) + (t0 - 
                                                  1)
        last <- max(dat$year[dat$rec == 1]) + (t0 - 1)
        firstMod <- t0
        lastMod <- tn
        yrs <- sort(unique(dat$year[dat$rec == 1]), decreasing = FALSE)
        gaps <- NULL
        yrs_data <- length(yrs) # number of years with data
        if (length(yrs) > 1) {
          for (i in (1:length(yrs) - 1)) {
            gaps <- c(gaps, yrs[i + 1] - yrs[i])
            
          }
        }
        if (!is.null(gaps)) {
          gap <- max(gaps)
        }
        else {
          gap <- 1
        }
        out2 <- data.frame(species, nRec_glob, nRec_reg, 
                           first, last, gap, firstMod, lastMod, yrs_data)
        rot <- attr(out_meta, "metadata")$analysis$spp_Metrics # LB removed 'as.data.frame', makes it an empty df instead of NULL
        # which causes next step to be missed
        if (is.null(rot)){
          rot <- data.frame(median = NA, P90 = NA, visits_median = NA, 
                            visits_P90 = NA, prop_list_one = NA, prop_repeats_grp = NA, 
                            prop_abs = NA)
        }
        #if(!is.null(rot)){ #LB
        rot <- as.data.frame(rot) #LB
        #} #LB
        rot$EqualWt <- ifelse(rot$prop_abs >= 0.99, rot$P90 >= 
                                3.1, rot$P90 >= 6.7)
        rot$HighSpec <- ifelse(rot$prop_abs >= 0.958, 
                               rot$P90 >= 9.5, rot$P90 >= 29)
        out2 <- cbind(out2, rot)
        print(paste("Sampled:", species))
        return(list(out1, out2))
      }
      else {
        print(paste("Error loading model:", species))
        return(NULL)
      }
    }
    else {
      if (!is.na(nRec) & !nRec >= minObs) 
        print(paste("Dropped (lack of observations):", 
                    species))
      else if (!is.na(nRec) & !region %in% c(out_meta$regions, 
                                             names(out_meta$region_aggs))) 
        print(paste("Dropped (region or region_aggs not present for species):", 
                    species))
      else print(paste("Error loading model:", species))
      return(NULL)
    }
    if(!is.null(out_meta) & (is.null(out_meta$region_aggs) | length(out_meta$region_aggs) == 0)){ # LB 
      REGION_IN_Q <- "psi.fs.r_UK"     # LB: for models with no region_agg data, rename as UK
      print("no region or region_aggs in model: assuming occupancy estimates are for UK")
    }#LB
  }
  if (parallel){ 
    outputs <- parallel::mclapply(keep, mc.cores = n.cores, 
                                  combineSamps, minObs = minObs)
  }
  else outputs <- lapply(keep, combineSamps, minObs = minObs)
  if (parallel) 
    samp_post <- parallel::mclapply(outputs, mc.cores = n.cores, 
                                    function(x) y <- x[[1]])
  else samp_post <- lapply(outputs, function(x) y <- x[[1]])
  samp_post <- do.call("rbind", samp_post)
  if (parallel) 
    meta <- parallel::mclapply(outputs, mc.cores = n.cores, 
                               function(x) y <- x[[2]])
  else meta <- lapply(outputs, function(x) y <- x[[2]])
  meta <- do.call("rbind", meta)
  meta <- data.frame(Species = meta$species, n_obs_global = meta$nRec_glob, 
                     n_obs_regional = meta$nRec_reg, min_year_data = meta$first, 
                     max_year_data = meta$last, min_year_model = meta$firstMod, 
                     max_year_model = meta$lastMod, gap_start = 0, gap_end = 0, 
                     gap_middle = meta$gap, yrs_data = meta$yrs_data, rot_median = meta$median, rot_P90 = meta$P90, 
                     rot_visits_median = meta$visits_median, rot_visits_P90 = meta$visits_P90, 
                     rot_prop_list_one = meta$prop_list_one, rot_prop_repeats_grp = meta$prop_repeats_grp, 
                     rot_prop_abs = meta$prop_abs, rot_EqualWt = meta$EqualWt, 
                     rot_HighSpec = meta$HighSpec)
 
  colnames(meta) <- paste0(colnames(meta), "_r_", gsub("psi.fs.r_", 
                                                       "", REGION_IN_Q))
  if (write == TRUE) {
    save(samp_post, file = paste(output_path, group_name, 
                                 "_all_spp_sample_", sample_n, "_post_", REGION_IN_Q, 
                                 ".rdata", sep = ""))
  }
  return(list(samp_post, meta))
}