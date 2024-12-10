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
library(BRCindicators) # BRCindicators: abundance and occupancy indicators
library(lubridate)
library('rslurm')
library(nimble) 
library(coda)
library(parallel)

# From DataLabs
# library(occNimble) # nicks nimble implementation

# Same parameters as Nick's run
maxSp <- 999 #5
yps <- 2 #4
n.iter <- 32e3 # 34e3
n.burn <- 30e3
min.Recs <- 50 # number of records per species for inclusion
inclStateRE <- TRUE
inclPhenology <- FALSE
ListLen <- "cat"
yps <- 2 #4
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

# Generate the job name with the current date
jobname <- paste0('dylcar_explore_occ_run_NIMBLE_', format(Sys.Date(), "%d_%m_%Y"))

# Detect the number of cores
cores_for_clustering <- 15 # We should explore here whether it is worthwile applying multiple species to each core if the species are processed easily within the wall time

# Calculate the total number of species
total_species <- length(formattedData$sp2incl)

# Calculate the number of nodes required
n_nodes_required <- ceiling(total_species / cores_for_clustering)

# Create a species-to-node mapping
species_to_node <- rep(1:n_nodes_required, each = cores_for_clustering, length.out = total_species)

# Define the Slurm function
slurm_NIMBLE_occ <- function(node_id) {

  transferable_outputs = initialise_NIMBLE()

    # Filter species for the given node
    species_indices <- which(species_to_node == node_id)
    sub_species_data <- transferable_outputs$data[[1]][species_indices, , drop = FALSE]

    # Process species in parallel within the node
    yearEff <- parallel::mclapply(species_indices, function(sp_i) {
        single_species_model(sp = sp_i,
                              spDat = list(sub_species_data[sp_i, , drop = FALSE]),
                              dataSumm = transferable_outputs$dataSumm,
                              n.iter = transferable_outputs$n.iter,
                              n.burn = transferable_outputs$n.burn,
                              n.thin = n.thin,
                              n.chain = n.chain,
                              Cmodel = transferable_outputs$Cmodel,
                              CoccMCMC = transferable_outputs$CoccMCMC,
                              mon2 = transferable_outputs$mon2)
    }, mc.cores = cores_for_clustering)

    names(yearEff) <- dimnames(transferable_outputs$dataSumm$occMatrix)[[1]][species_indices]
    attr(yearEff, "modelCode") <- transferable_outputs$model$getCode()

    saveRDS(yearEff, file = paste0("yearEff_node_", node_id, ".rds"))
}

# Create the Slurm job
sjob <- slurm_apply(
    f = slurm_NIMBLE_occ,
    params = data.frame(node_id = 1:n_nodes_required),
    jobname = jobname,
    nodes = n_nodes_required, # ask for extra core
    cpus_per_node = cores_for_clustering,
    submit = TRUE,
    global_objects = c("initialise_NIMBLE", "defineModel_SS", "species_to_node", "n_nodes_required", 
                       "cores_for_clustering", "n.thin", "n.chain", "single_species_model", "n_nodes_required", "total_species", "formattedData", "dataSumm", "maxSp", "n.iter", "n.burn", "inclStateRE", "inclPhenology", "ListLen", "all.Pars", "ListLen"),
    slurm_options = list(time = '10:00:00', mem = 20000, partition = 'par-single', error = '%a.err')
)
