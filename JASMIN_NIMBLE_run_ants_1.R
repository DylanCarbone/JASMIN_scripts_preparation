# setup the session and R ...
module load jaspy # For loading python workspace ??
module load jasr # For loading R workspace ??
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

sjob <- slurm_call(f = initialise_NIMBLE,
                   jobname = jobname,
                   global_objects = c("formattedData", "dataSumm", "inclStateRE", "ListLen", "inclPhenology", "parallel", "all.Pars", "n.iter", "n.burn", "maxSp", "defineModel_SS"),
                   slurm_options = list(time = '10:00:00',
                                        mem = 20000,
                                        partition = 'short-serial',
                                        error = '%a.err'))

# For species model run
  # n.thin = 5,
  # n.chain = 3,
