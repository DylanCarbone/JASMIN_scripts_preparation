# Reformat the occurrence data into a format suitable for submitting to occDetFunc

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

# Same parameters as Nick's run
maxSp <- 999 #5
yps <- 2 #4
n.iter <- 32e3 # 34e3
n.burn <- 30e3
saveData <- TRUE
min.Recs <- 50 # number of records per species for inclusion
# NOTE: I have set the maximum tasks to 7 in the rslurm templates for 7 cores


# in case we need to run scripts in user directory
setwd("dylcar_explore_occ_user")

# Source modified functions
for (file in list.files("modified_functions", pattern = "*.R")) {
  source(file.path("modified_functions", file))
}

# Read in the occurrence data
ants_data <- readRDS("monad_occupancy_dataset_ants.rds")

ants_visit_data <- formatOccData(taxa = ants_data$tik,
                                 site = ants_data$GRIDREF,
                                 survey = ants_data$lower)

saveRDS(ants_visit_data, "ants_visit_data.rds")

# Load formatted taxon group data
ants_visit_data <- readRDS(file = "ants_visit_data.rds")

# Load region data - This was used in Kath's scripts but is excluded for the purpose of this workflow
# reg_data <- read.csv("UK_grid_refs_monads.csv", stringsAsFactors = FALSE)
# reg_data <- reg_data[, 1:5]
# names(reg_data)[1] <- "site"
# region_aggs <- list(GB = c('ENGLAND', 'WALES', 'SCOTLAND'), UK = c('ENGLAND', 'WALES', 'SCOTLAND', 'NORTHERN_IRELAND'))

# Initialize the log CSV file
log_df <- data.frame(
  taxa = character(),
  run_type = character(),
  parallel = logical(),
  start_time = character(),
  taxa_name = character(),
  node_start_time = character(),
  node_end_time = character(),
  node_run_time_hours = numeric(),
  output_object_size_mb = numeric(),
  message = character(),
  stringsAsFactors = FALSE
)

# Save the initial empty log file
write.csv(log_df, "runtime_log.csv", row.names = FALSE)

# Create slurm function with logging
slurm_occDetFunc <- function(taxa_name) {
  node_start_time <- Sys.time()
  
  log_entry <- data.frame(
    taxa = "ant", # Adding the fixed column 'taxa' with the value "ant"
    run_type = "JASMIN",
    parallel = TRUE,
    start_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    taxa_name = taxa_name,
    node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
    node_end_time = NA,
    node_run_time_hours = NA,
    output_object_size_mb = NA,
    message = NA,
    stringsAsFactors = FALSE
  )
  
  tryCatch({

    out <-  occDetFunc(taxa_name = as.character(taxa_name),
                output_dir = "occ_outputs",
                occDetdata = ants_visit_data$occDetdata,
                spp_vis = ants_visit_data$spp_vis,
                write_results = saveData,
                n_chains = 3,
                n_iterations = n.iter,
                burnin = n.burn,
                thinning = 6,
                nyr = yps,
                criterion = min.Recs,
                modeltype = c('ranwalk', 'halfcauchy', 'catlistlength'),
                return_data = FALSE,
                provenance = "Data generated in a trial of sparta performance and BRC workflow performance in JASMIN. PLEASE DO NOT SHARE THE OUPUTS OF THIS WORKFLOW")
    
    log_entry$output_object_size_mb = as.numeric(object.size(out)/(1024^2))
    
    node_end_time <- Sys.time()
    node_run_time <- node_end_time - node_start_time
    
    # Log successful completion
    log_entry$node_end_time <- format(node_end_time, "%Y-%m-%d %H:%M:%S")
    log_entry$node_run_time_hours <- as.numeric(node_run_time)
    log_entry$message <- "Success"
    
  }, error = function(e) {
    # Log error
    log_entry$message <- as.character(e$message)
  })
  
  # Append the result to the log CSV file, ensuring column titles are not repeated
  write.table(log_entry, "../runtime_log.csv", sep = ",", row.names = FALSE, 
              col.names = !file.exists("../runtime_log.csv"), append = TRUE)
  
  return(NULL)
}

# Create roster of parameters
pars <- data.frame(taxa_name = as.character(names(ants_visit_data[['spp_vis']])[-1]))

#Adjust length by number of species
pars <- pars %>% slice(1:ifelse(maxSp > nrow(pars), nrow(pars), maxSp))

# Generate the job name with the current date
jobname <- paste0('dylcar_explore_occ_run_', format(Sys.Date(), "%d_%m_%Y"))

# Create the job using rslurm
sjob <- slurm_apply(f = slurm_occDetFunc,
                    params = pars,
                    jobname = jobname,
                    nodes = nrow(pars),
                    cpus_per_node = 1,
                    submit = TRUE,
                    global_objects = c('ants_visit_data', "saveData", "n.iter", "n.burn", "yps", "min.Recs"),
                    slurm_options = list(time = '10:00:00', 
                                         mem = 20000,
                                         partition = 'short-serial',
                                         error = '%a.err'))

