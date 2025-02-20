# setup the session and R ...
module load jaspy # For loading python workspace
module load jasr # For loading R workspace
R

# remotes::install_github("BiologicalRecordsCentre/sparta")

# remotes::install_github("biologicalrecordscentre/BRCindicators")

# remotes::install_github("BiologicalRecordsCentre/wrappeR")

library(remotes) # remotes: package installation
library(withr) # withr: package installation
library(plyr) # plyr: data manipulation
library(dplyr) # dplyr: data manipulation
library(tidyr) # tidyr: data manipulation
library(reshape2) # reshape2: data manipulation
library(rjags) # rjags: occupancy models
library(R2jags) # R2jags: occupancy models
library(jagsUI)
library(ggplot2) # ggplot2: plotting
library(gridExtra) # gridExtra: plotting
library(cowplot) # cowplot: plotting
library(readr)
library(stringr)
library(sparta) # sparta: occupancy models
library(lubridate)
library('rslurm')
library(parallel)
library(nimble)

# consider memory weighed against time constraints
cores_for_clustering <- 6 # We should explore here whether it is worthwile applying multiple species to each core if the species are processed easily within the wall time
species_to_core <- 2

# Generate the job name with the current date
jobname <- paste0('dylcar_explore_occ_run_NIMBLE_6_c_2_spc', "_", format(Sys.Date(), "%d_%m_%Y"))

# Same parameters as Nick's run
maxSp <- 999 #5
yps <- 2 #4
n.iter <- 32e3 # 32e3
n.burn <- 30e3 # 30e3
min.Recs <- 50 # number of records per species for inclusion
inclStateRE <- TRUE
inclPhenology <- FALSE
ListLen <- "cat"
all.Pars <- c("sd.psi", "lam.0", "mu.alpha", "sd.alpha", "gamma.1", "gamma.2")#, "psi.fs")
n.thin = 5
n.chain = 3

# NOTE: I have set the maximum tasks to 7 in the rslurm templates for 7 cores

# in case we need to run scripts in user directory
setwd("dylcar_explore_occ_user")

# Source modified functions
for (file in list.files("nimble_implementation_functions", pattern = "*.R")) {
  source(file.path("nimble_implementation_functions", file))
}

# read in the ant data
ants_data <- readRDS("monad_occupancy_dataset_ants.rds")

# convert the names to our desired format
inData <- ants_data %>%
  rename(species = tik,
        siteID = GRIDREF,
        survey = lower_date) %>%
  mutate(year = as.numeric(format(as.Date(survey, format="%d/%m/%Y"),"%Y")))

# prepare the data for the model (includes removing species found on few sites)
formattedData <- formatData(inData, 
                            inclPhenology = inclPhenology,
                            ListLen = ListLen,
                            minYrPerSite = yps,
                            minSite = 1, # sites with records per species
                            minRecs = min.Recs) # records per species

# Summarise the data - need the Nimble formatted data for this
dataSumm <- with(formatData(inData, 
                            verbose = FALSE,
                            format = "Nimble"), 
                  summariseData(obsData, dataConstants))

# Calculate the total number of species
total_species <- length(formattedData$sp2incl)

# Calculate the number of nodes required
n_nodes_required <- ceiling(total_species / (cores_for_clustering * species_to_core))

# Create a species-to-node mapping
species_to_node <- rep(1:n_nodes_required, each = cores_for_clustering, length.out = total_species)

# Define the Slurm function
slurm_NIMBLE_occ <- function(node_i) {

  node_start_time <- Sys.time()

  transferable_outputs <- initialise_NIMBLE()

  initialisation_run_time <- as.numeric((Sys.time() - node_start_time) / 60)

  # Filter species for the given node
  species_indices <- which(species_to_node == node_i)

  # Group species into chunks based on species_to_core
  species_groups <- split(species_indices, ceiling(seq_along(species_indices) / species_to_core))

  logs <- parallel::mclapply(species_groups, function(group) {
    core_start_time <- Sys.time()

    log_entries <- list()

    for (core_i in group) {
      species_name <- transferable_outputs$dataSumm$stats$species[core_i]

      # Initialise the log entry
      log_entry <- data.frame(
        taxa_group = "Ants",
        species_name = species_name,
        JASMIN = TRUE,
        queue = "par-single",
        n_nodes_requested = n_nodes_required,
        n_node_cores = cores_for_clustering,
        NIMBLE_initialisation_run_time = initialisation_run_time,
        node_iteration = node_i,
        node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
        node_end_time = NA,
        node_run_time = NA,
        core_iteration = core_i,
        core_start_time = format(core_start_time, "%Y-%m-%d %H:%M:%S"),
        core_end_time = NA,
        core_run_time = NA,
        output_object_size_GB = NA,
        stringsAsFactors = FALSE
      )

      tryCatch({
        # Run the single species model
        sp_output <- single_species_model(
          sp = core_i,
          spDat = list(transferable_outputs$data[[1]][rownames(transferable_outputs$data[[1]]) == species_name, , drop = FALSE]),
          dataSumm = transferable_outputs$dataSumm,
          n.iter = transferable_outputs$n.iter,
          n.burn = transferable_outputs$n.burn,
          n.thin = n.thin,
          n.chain = n.chain,
          Cmodel = transferable_outputs$Cmodel,
          CoccMCMC = transferable_outputs$CoccMCMC,
          mon2 = transferable_outputs$mon2
        )

        # Save the species output
        saveRDS(sp_output, paste0(species_name, ".rds"))
        log_entry$output_object_size_GB <- as.numeric(object.size(sp_output)) / (1024^3)
      }, error = function(e) {
        # Save error message to a text file
        error_message <- paste0("Error for species: ", species_name, "\n", conditionMessage(e))
        writeLines(error_message, paste0(species_name, "_error.txt"))
      })

      core_end_time <- Sys.time()

      log_entry$core_end_time <- format(core_end_time, "%Y-%m-%d %H:%M:%S")
      log_entry$core_run_time <- as.numeric(difftime(core_end_time, core_start_time, units = "hours"))

      log_entries[[length(log_entries) + 1]] <- log_entry

      # Clean up memory for the current species
      gc()
    }

    do.call(rbind, log_entries)
  }, mc.cores = cores_for_clustering)

  node_log <- do.call(rbind, logs)

  node_end_time <- Sys.time()

  node_log$node_end_time <- format(node_end_time, "%Y-%m-%d %H:%M:%S")
  node_log$node_run_time <- as.numeric(difftime(node_end_time, node_start_time, units = "hours"))

  write.csv(node_log, paste0("node_", node_i, "_log.csv"))
}

# Create the Slurm job
sjob <- slurm_apply(
    f = slurm_NIMBLE_occ,
    params = data.frame(node_i = 1:n_nodes_required),
    jobname = jobname,
    nodes = n_nodes_required + 1, # ask for extra core
    cpus_per_node = cores_for_clustering,
    submit = TRUE,
    global_objects = c("initialise_NIMBLE", "defineModel_SS", "species_to_node", "n_nodes_required", 
                       "cores_for_clustering", "n.thin", "n.chain", "single_species_model", "n_nodes_required", "total_species", "formattedData", "dataSumm", "maxSp", "n.iter", "n.burn", "inclStateRE", "inclPhenology", "ListLen", "all.Pars", "species_to_core"),
    slurm_options = list(time = '24:00:00', mem = 40000, partition = 'par-single', error = '%a.err')
)
