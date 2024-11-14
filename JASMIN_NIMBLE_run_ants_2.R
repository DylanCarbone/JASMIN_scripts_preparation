# setup the session and R ...
module load jaspy # For loading python workspace ??
module load jasr # For loading R workspace ??

# Check your job has completed
squeue --user dylcar

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

# From DataLabs
# library(occNimble) # nicks nimble implementation

# Same parameters as Nick's run
yps <- 2 #4
n.iter <- 32e3 # 34e3
n.burn <- 30e3
min.Recs <- 50 # number of records per species for inclusion
inclPhenology <- FALSE
ListLen <- "cat"
yps <- 2 #4
n.thin = 5
n.chain = 3

# NOTE: I have set these maximum tasks to 7 in the rslurm templates for 7 cores

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

jobname <- "_rslurm_dylcar_explore_occ_run_NIMBLE_11_11_2024"
#jobname <- paste0('dylcar_explore_occ_run_NIMBLE_', format(Sys.Date(), "%d_%m_%Y")) # may have to change if ran on a different day

# May take a couple minutes to load...
# read in the transferable_outputs
transferable_outputs <- readRDS("transferable_outputs.rds")
transferable_outputs <- readRDS(file.path(jobname, "transferable_outputs.rds"))

# Run for a single species
species_run <- function(sp_i){

  yearEff <- single_species_model(sp = sp_i,
                        spDat = lapply(transferable_outputs$obsData, function(x) x[sp_i,]),
                        dataSumm = transferable_outputs$dataSumm,
                        n.iter = transferable_outputs$n.iter,
                        n.burn = transferable_outputs$n.burn,
                        n.thin = n.thin,
                        n.chain = n.chain,
                        transferable_outputs$Cmodel, transferable_outputs$CoccMCMC,
                        mon2 = transferable_outputs$mon2)

  print("function runs")

  saveRDS(yearEff, file = paste0(dimnames(formattedData$obsData$y)[[1]][1:transferable_outputs$nSpMod], "_model_run.rds"))

  # attr(yearEff, "modelCode") <- model$getCode()  # Should append this to the final outputs
}

test_sp =  data.frame(sp_i = 1:transferable_outputs$nSpMod)$sp_i[1]

species_run(test_sp)

# # Create the job using rslurm that compiles the initial model
# sjob <- slurm_apply(f = species_run,
#                     params = data.frame(sp_i = 1:transferable_outputs$nSpMod),
#                     jobname = jobname,
#                     nodes = transferable_outputs$nSpMod, # Assuming this cannot be parallelised
#                     cpus_per_node = 1, # We need to look into whether this improves run time
#                     submit = TRUE,
#                     global_objects = c("formattedData", "transferable_outputs", "n.thin", "n.chain"),
#                     slurm_options = list(time = '10:00:00', 
#                                          mem = 20000,
#                                          partition = 'short-serial',
#                                          error = '%a.err'))
