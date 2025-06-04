library(base, quietly = TRUE)
library(methods, quietly = TRUE)
library(datasets, quietly = TRUE)
library(utils, quietly = TRUE)
library(grDevices, quietly = TRUE)
library(graphics, quietly = TRUE)
library(stats, quietly = TRUE)
library(parallel, quietly = TRUE)
library(rslurm, quietly = TRUE)
library(nimble, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tibble, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(readr, quietly = TRUE)
library(purrr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(forcats, quietly = TRUE)
library(lubridate, quietly = TRUE)
library(MCMCvis, quietly = TRUE)
library(Matrix, quietly = TRUE)

rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

allSpecies = readRDS("allSpecies.RDS")
nimbleConstants = readRDS("nimbleConstants.RDS")
mynimbleCode = readRDS("mynimbleCode.RDS")
visit_df = readRDS("visit_df.RDS")
ni = readRDS("ni.RDS")
nb = readRDS("nb.RDS")
nt = readRDS("nt.RDS")
n.chains = readRDS("n.chains.RDS")
species_counts = readRDS("species_counts.RDS")

NIMBLE_run_parallel = function(species_i){

  print(paste("species_i:", species_i))
  
  node_start_time = Sys.time()

  myspecies <- allSpecies[species_i]
  print(myspecies)
  
  # Match visits and replace NA with 0
  visit_df$Species <- ifelse(visit_df$visit %in% names(species_counts), species_counts[visit_df$visit], 0)

  nimbleData <- list(y = as.numeric(visit_df$Species))
  
  #specify inits
  z_inits <- visit_df %>%
    group_by(siteIndex, yearIndex) %>%
    summarise(Species = max(Species, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = yearIndex, values_from = Species, values_fill = 0) %>%
    select(-siteIndex) %>%
    as.matrix()

  #ignore warning message
  z_inits[is.na(z_inits)] <- 0
  z_inits[!z_inits %in% c(0,1)]<- 0

  nimbleInits <- list(z = z_inits)

  # Define log file path
  log_file <- paste0("nimble_run_log.txt")
  
  myModel <- nimbleModel(code = mynimbleCode,
                          data = nimbleData,
                          constants = nimbleConstants,
                          inits = nimbleInits)

  write("Model initialized", file = log_file, append = TRUE)

  CmyModel <- compileNimble(myModel)

  write("Model compiled", file = log_file, append = TRUE)

  myMCMC <- buildMCMC(CmyModel, monitors = c('sd.s','year.fixed','prop.p','beta.p'))

  write("MCMC built", file = log_file, append = TRUE)

  CmyMCMC <- compileNimble(myMCMC)

  write("MCMC compiled", file = log_file, append = TRUE)

  initialisation_run_time <- as.numeric(difftime(Sys.time(), node_start_time, units = "hours"))
  
  parallel::mclapply(1:n.chains, function(chain){

    chain_log_file <- paste0("chain_", chain, "_MCMC_log.txt")

    results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1)

    write("MCMC run completed", file = chain_log_file, append = TRUE)

    #save output
    temp <- MCMCsummary(object = results, round = 5)
    saveRDS(temp, file=paste0("mcmc.summary_",myspecies,".rds"))

    return(NULL)
  })
  
  log_entry <- data.frame(
    taxa_group = "Ants with 3 parallel threads",
    species_name = myspecies,
    JASMIN = TRUE,
    queue = "par-single",
    n_nodes_requested = length(allSpecies),
    NIMBLE_initialisation_run_time = initialisation_run_time,
    node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
    node_end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    node_run_time = as.numeric(difftime(Sys.time(), node_start_time, units = "hours")),
    stringsAsFactors = FALSE
  )
  write.csv(log_entry, paste0(myspecies, "_log.csv"))
}

NIMBLE_run_parallel(rslurm_id + 1)