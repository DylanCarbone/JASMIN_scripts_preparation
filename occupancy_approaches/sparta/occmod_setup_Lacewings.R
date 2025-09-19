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
setwd("/gws/ssde/j25b/ceh_generic/dylcar/indicators_2025/Lacewings_Sp")

# Format Lacewings data - remove already completed species
data = read.csv("Data_LAC.csv") %>%
rename(date = Start_Date, gridref = Source_Location, species = Taxon_out) %>%
mutate(species = gsub("/", "_", species), 
date = ymd(date), year = year(date), gridref)

# Use same name as previos runs/specify name of new run
jobname = "_rslurm_Lacewings"

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
                                
saveRDS(visitData, "Lacewings_Sp_filtered_start_years.rds")

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

# NB: This deletes old logs inside the job folder
old_job_err_files <- list.files(jobname, pattern = "\\.err$", full.names = TRUE)
old_job_out_files <- list.files(jobname, pattern = "\\.out$", full.names = TRUE)

# delete the files
file.remove(old_job_err_files)
file.remove(old_job_out_files)

# log how many records there are for each years
if(file.exists("species_number_tally.csv")){
  file.remove("species_number_tally.csv")
}

species_number_tally = data %>% filter(species %in% species_to_run) %>%
group_by(species, year) %>%
summarise(n_entries = n())

# Just to stop tibble concatenating it
as.data.frame(species_number_tally)

write.csv(species_number_tally, "species_number_tally.csv")


# Create roster
pars <- data.frame(taxa_name = species_to_run)

sjob <- slurm_apply(f = slurm_occDetFunc,
                    params = pars, 
                    # to change for each run
                    jobname = 'Lacewings',   #!!!
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
