# setup the session and R ...
module load jasr

R

# Set working directory
setwd("dylcar_explore_occ_user")

# load libraries
require('rslurm')
require('sparta')
require('reshape2')
require(dplyr)
require(rjags)
require(R2jags)
require(jagsUI)
require('rslurm')
require('sparta')
require('reshape2')
require(ggplot2)

# Load data --------------------------------------------------------------------

# species occurrence records - essentially a list of species / date / 1km grids
visitData <- readRDS(file = "monad_occupancy_dataset_ants.rds")

taxa_group = "Ants"

# spatial information ----------------------------------------------------------

#sometimes we run model that include a factor for region/country
# load region data
reg_data <- read.csv("UK_grid_refs_monads.csv", stringsAsFactors = FALSE) # From the folder I have saved the GB_grid_refs.csv file in the MobaXterm session folders
reg_data <- reg_data[,1:5]
names(reg_data)[1] <- "site"
region_aggs <- list(GB = c('ENGLAND','WALES','SCOTLAND'), UK = c('ENGLAND','WALES','SCOTLAND','NORTHERN_IRELAND'))

# Function ---------------------------------------------------------------------

slurm_occDetFunc <- function(species_i){

    node_start_time = Sys.time()

    taxa_name = unique(visitData$tik)[species_i]

    formatted_data <- formatOccData(taxa = visitData$tik,
                             site = visitData$GRIDREF,
                             survey = visitData$lower_date)

    if (!dir.exists("results")){
      dir.create("results")
    }

    out <- occDetFunc(taxa_name = as.character(taxa_name),
                      occDetdata = formatted_data$occDetdata,
                      spp_vis = formatted_data$spp_vis,
                      write_results = TRUE,
                      output_dir = "results",
                      n_chains = 3,
                      n_iterations = 32000,
                      burnin = 30000,
                      thinning = 6,
                      nyr = 2,
                      modeltype = 'sparta',
                      #regional_codes = reg_data,
                      #region_aggs = region_aggs,
                      return_data = FALSE,
                      seed = 123,
                      additional.parameters = "a",
                      allowSitesMultiRegions = TRUE,
                      rem_aggs_with_missing_regions=FALSE,
                      provenance = "test")

        # Extract trend data for a species
    plot = plot(out)

    log_entry <- data.frame(
    taxa_group = "Ants", ### HERE
    species_name = taxa_name,
    JASMIN = TRUE,
    queue = "long-serial",
    n_nodes_requested = length(unique(visitData$tik)),
    node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
    node_end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    node_run_time = as.numeric(difftime(Sys.time(), node_start_time, units = "hours")),
    n_start_val = NA,
    stringsAsFactors = FALSE
  )

    if (!dir.exists("logs")){
      dir.create("logs")
    }

    write.csv(log_entry, paste0("logs/", taxa_name,  "_log.csv"))
    
    if (!dir.exists("plots")){
      dir.create("plots")
    }
    
    ggsave(paste0("plots/", taxa_name, ".png"), plot = plot)

    return(NULL)
}

# Create roster ----------------------------------------------------------------

# define which species we are going to run the model for 

pars <- data.frame(species_i = 1:length(unique(visitData$tik)))

# create slurm job --------------------------------------------------------------

# Create the job scipt and the R script needed to run the process on 
# lotus using slurm. Note: you can edit the templates used. These are
# found in the slurm folder in your R library (run '.Library' to find).
# You will need to add the command to load jaspy: module add jaspy

# Generate the job name with the current date
jobname <- paste0('dylcar_explore_occ_run_SPARTA_', toupper(taxa_group), "_", format(Sys.Date(), "%d_%m_%Y"))
dir.create(paste0("_rslurm_", jobname))

# Create the job
sjob <- slurm_apply(f = slurm_occDetFunc,
                    params = pars,
                    jobname = jobname, # modify for group e.g., Ladybirds_2022
                    nodes = nrow(pars),
                    cpus_per_node = 1,
                    submit = TRUE,
                    global_objects = c('visitData', 'reg_data', 'region_aggs'),
  slurm_options = list(time = "24:00:00", mem = 30000, error = "%a.err",
  account = "ceh_generic", partition = "standard", qos = "long"))

# end --------------------------------------------------------------------------
