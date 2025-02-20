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
occMatrix = readRDS("occMatrix.RDS")
ni = readRDS("ni.RDS")
nb = readRDS("nb.RDS")
nt = readRDS("nt.RDS")
n.chains = readRDS("n.chains.RDS")
species_to_node = readRDS("species_to_node.RDS")

NIMBLE_run_parallel = function(node_i){

    node_start_time <- Sys.time()

    node_species <- allSpecies[which(species_to_node == node_i)]

    #identify the first node species
    first_species = node_species[1]

    visit_df$first_species <- occMatrix[,first_species]
    nimbleData <- list(y = as.numeric(visit_df$first_species))

    z_inits <- reshape2::acast(visit_df, siteIndex ~ yearIndex, value.var="first_species", fun = max)
    z_inits[is.na(z_inits)] <- 0
    z_inits[!z_inits %in% c(0,1)]<- 0
    
    nimbleInits <- list(z = z_inits)

    # Define log file path
    log_file <- paste0("nimble_initialisation_run_log.txt")
  
    # 
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

    initialisation_run_time <- as.numeric((Sys.time() - node_start_time) / 60)

    parallel::mclapply(1:length(node_species), function(core_i){

        core_species = node_species[core_i]

        # Define log file path
        species_log_file <- paste0(core_species, "_model_fitting_log.txt")

        write(paste0("handling: ", core_species), file = species_log_file, append = TRUE)

        visit_df$Species <- occMatrix[, core_species]

        nimbleData <- list(y = as.numeric(visit_df$Species))
  
        z_inits <- reshape2::acast(visit_df, siteIndex ~ yearIndex, value.var="Species", fun = max)
        z_inits[is.na(z_inits)] <- 0
        z_inits[!z_inits %in% c(0,1)]<- 0
        
        nimbleInits <- list(z = z_inits)

        CmyModel$setData(nimbleData)
        write("Data set successfully", file = species_log_file, append = TRUE)

        CmyModel$setInits(list(z = z_inits))
        write("Initials set successfully", file = species_log_file, append = TRUE)

        results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = 1)
        write("MCMC run completed", file = species_log_file, append = TRUE)

        #save output
        temp <- MCMCsummary(object = results, round = 5)
        saveRDS(temp, file=paste0("mcmc.summary_",core_species,".rds"))

        return(NULL)
    })

    log_entry <- data.frame(
        taxa_group = "Ants with multiple species per core",
        species_name = paste(node_species, collapse = " and "),
        JASMIN = TRUE,
        queue = "par-single",
        n_nodes_requested = length(unique(species_to_node)),
        NIMBLE_initialisation_run_time = initialisation_run_time,
        node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
        node_end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        node_run_time = as.numeric(difftime(Sys.time(), node_start_time, units = "hours")),
        stringsAsFactors = FALSE
    )
  write.csv(log_entry, paste0(node_i, "_log.csv"))
}

NIMBLE_run_parallel(rslurm_id + 1)
