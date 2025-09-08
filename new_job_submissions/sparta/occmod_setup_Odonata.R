## setting up the slurm scripts to run models on jasmin
## run on NERC supercomputer jasmin
## Author: Francesca Mancini

## Bees ----

rm(list = ls())

# load libraries
require('BRCmap')
require('rslurm')
require('sparta')
require('reshape2')
require('lubridate')
require('dplyr')

# Load data
# to change for each run # !!!
setwd("/gws/nopw/j04/ceh_generic/dylcar/indicators_2025//Odonata_Sp")

# Format hoverfly data - remove already completed species
data = read.csv("odonata_cleaned_data_2025.csv") %>% mutate(species = gsub("/", "_", species), 
date = ymd(date), year = year(date))

# Use same name as previos runs/specify name of new run
jobname = "_rslurm_Dragonfly"

# Filter out species that have already been modelled
skip_sp = gsub(".rdata", "", basename(list.files(jobname, pattern = "*.rdata", recursive = TRUE)))

# Filter sites with >= 2 unique years before computing first_complete_year
data_temp <- data %>%
  group_by(gridref) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

first_complete_year <- {
  y <- max(data_temp$year)
  while ((y - 1) %in% data_temp$year) y <- y - 1
  y
}

if(first_complete_year != min(data$year)){ paste("first year is now:", first_complete_year)}

data <- data %>% filter(year >= first_complete_year)

# run the model with these data for one species
visitData <- formatOccData(taxa = data$species,
                                site = data$gridref,
                                survey = data$date)
                                
saveRDS(visitData, "odonata_Sp_filtered_start_years.rds")

reg_data <- readRDS("../sq1km_UK_regions.rds")

region_aggs <- list(UK = as.character(c("ENGLAND","SCOTLAND", 
                                        "WALES", "NORTHERN_IRELAND")),
                    GB = as.character(c("ENGLAND","SCOTLAND", "WALES")))

# Define function that loops through species
slurm_occDetFunc <- function(taxa_name){
  
  out <- sparta::occDetFunc(taxa_name = as.character(taxa_name),
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
                    provenance = "2025 indicator update")
  return(NULL)
}

species_to_run = data %>% filter(!species %in% skip_sp) %>% pull(species) %>% unique()

# Create roster
pars <- data.frame(taxa_name = species_to_run)

sjob <- slurm_apply(f = slurm_occDetFunc,
                    params = pars, 
                    # to change for each run
                    jobname = 'Hoverflies',   #!!!
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
