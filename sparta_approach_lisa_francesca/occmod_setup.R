## setting up the slurm scripts to run models on jasmin
## run on NERC supercomputer jasmin
## Author: Francesca Mancini

rm(list = ls())

# load libraries
require('rslurm')
require('sparta')
require('reshape2')

# Load data
# to change for each run
setwd("/gws/nopw/j04/ceh_generic/framan/indicators_2025/wasps")

# to change for each run
visitData <- readRDS(file = './Was_Sp.rds')

reg_data <- readRDS("../sq1km_UK_regions.rds")


region_aggs <- list(UK = as.character(c("ENGLAND","SCOTLAND", 
                                        "WALES", "NORTHERN_IRELAND")),
                    GB = as.character(c("ENGLAND","SCOTLAND", "WALES")))

# Define function that loops through species
slurm_occDetFunc <- function(taxa_name){
  
  out <- occDetFunc(taxa_name = as.character(taxa_name),
                    occDetdata = visitData$occDetdata,
                    spp_vis = visitData$spp_vis,
                    write_results = TRUE,
                    n_chains = 3,
                    n_iterations = 32000,
                    burnin = 30000,
                    thinning = 6,
                    nyr = 2,
                    modeltype = c('ranwalk', 'halfcauchy', 'catlistlength'),
                    regional_codes = reg_data,
                    region_aggs = region_aggs,
                    return_data = FALSE,
                    seed = 123,
                    additional.parameters = "a",
                    allowSitesMultiRegions = FALSE,
                    rem_aggs_with_missing_regions=FALSE,
                    # to change for each run
                    provenance = "BWARS - 2025 indicator update") 
  return(NULL)
}

# Create roster
pars <- data.frame(taxa_name = as.character(names(visitData[['spp_vis']])[-1]))


sjob <- slurm_apply(f = slurm_occDetFunc,
                    params = pars, 
                    # to change for each run
                    jobname = 'wasps',
                    nodes = nrow(pars), 
                    cpus_per_node = 1, 
                    submit = TRUE,
                    global_objects = c('visitData', 'reg_data', 'region_aggs'),
                    slurm_options = list(account = 'ceh_generic',
                                         time = '119:59:00', 
                                         mem = 20000,
                                         partition = 'standard',
                                         qos = 'long',
                                         error = '%a.err'))
