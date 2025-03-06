# setup the session and R ...
module load jaspy # For loading python workspace ??
module load jasr # For loading R workspace ??

R

library(parallel)
library(rslurm)
library(nimble)
library(tidyverse)
library(MCMCvis)
library(dplyr)
library(Matrix)
library(rnrfa)

library(data.table)
library(unmarked)
library(ggplot2)
library(igr)

# Model parameters
min.Recs <- 50 # 50 # number of records per species for inclusion

# Set working directory
setwd("dylcar_explore_occ_user")

for (file in list.files("occti_functions", full.names = TRUE)){
source(file)
}

# Read and preprocess data
data <- readRDS("monad_occupancy_dataset_ants.rds") %>% 
  rename(Species = tik, SiteID = GRIDREF, Date = lower_date) %>%
  mutate(Date = as.Date(Date),
         yday = lubridate::yday(Date),
         Year = lubridate::year(Date))

# Rename columns to match expected function inputs
setnames(data, "SiteID", "Gridref")  # Rename SiteID to Gridref

#############################

gb_ir <- split(data, igr_is_valid(data$Gridref))
names(gb_ir) = ifelse(length(gb_ir) ==2, c("gb", "ir"), "gb")

# For now, I will use only british data
data = gb_ir$gb

############################

# Convert Date column and extract Year, Month, Week
data <- add_dates(data)

# Calculate list length (number of species observed per site-date)
data <- add_listL(data)

# Summarise species
speciesSummary <- data %>%
  group_by(Species) %>%
  summarise(nuSiteID = length(unique(Gridref)),
            nuRecs = length(Species)) %>%
  arrange(desc(nuRecs))

# Define species list
allSpecies <- sort(speciesSummary$Species[speciesSummary$nuRecs > min.Recs])

combinations <- c(
  "SA", "TA", "NA", "HA", "OA", "IA", "SB", "TB", "NB", "HB", "OB", "IB", "SC", "TC", "NC",
  "HC", "OC", "IC", "SD", "TD", "ND", "HD", "OD", "ID", "SE", "TE", "NE", "HE", "OE", "IE",
  "SF", "TF", "NF", "HF", "OF", "IF", "SG", "TG", "NG", "HG", "OG", "IG", "SH", "TH", "NH",
  "HH", "OH", "IH", "SJ", "TJ", "NJ", "HJ", "OJ", "IJ", "SK", "TK", "NK", "HK", "OK", "IK",
  "SL", "TL", "NL", "HL", "OL", "IL", "SM", "TM", "NM", "HM", "OM", "IM", "SN", "TN", "NN",
  "HN", "ON", "IN", "SO", "TO", "NO", "HO", "OO", "IO", "SP", "TP", "NP", "HP", "OP", "IP",
  "SQ", "TQ", "NQ", "HQ", "OQ", "IQ", "SR", "TR", "NR", "HR", "OR", "IR", "SS", "TS", "NS",
  "HS", "OS", "IS", "ST", "TT", "NT", "HT", "OT", "IT", "SU", "TU", "NU", "HU", "OU", "IU",
  "SV", "TV", "NV", "HV", "OV", "IV", "SW", "TW", "NW", "HW", "OW", "IW", "SX", "TX", "NX",
  "HX", "OX", "IX", "SY", "TY", "NY", "HY", "OY", "IY", "SZ", "TZ", "NZ", "HZ", "OZ", "IZ"
)

data <- data[data$Gridref %in% grep(paste0("^(", paste(combinations, collapse = "|"), ")"), data$Gridref, value = TRUE), ]

positions = rnrfa::osg_parse(grid_refs = data$Gridref, coord_system = "BNG")

data$North = positions$northing
data$East = positions$easting

data_range = range(data$Year)

occti_run = function(species_i){

    node_start_time = Sys.time()

    myspecies <- allSpecies[species_i]
    print(myspecies)

    # Run the single-species occupancy model
    occupancy_result <- fit_occ(
    spp = myspecies,
    obdata = data,
    occformula = "~ North + I(North^2) + East + I(East^2)",  # Quadratic spatial covariates
    detformula = "~ logLL + SEAS",  # Detection probability model
    covnames = c("East", "North"),  # Covariates used
    minyear = data_range[1],  # Start year
    maxyear = data_range[2],  # End year
    trendyears = data_range,  # Years for trend estimation
    nstart = 1,  # Number of initial values to test
    outputdir = "",  # Directory to save results
    engine = "C"
    )

    # Save the results
    saveRDS(occupancy_result, paste0(myspecies, "_occupancy_output.rds"))

    log_entry <- data.frame(
    taxa_group = "Butterflies", ### HERE
    species_name = myspecies,
    JASMIN = TRUE,
    queue = "long-serial",
    n_nodes_requested = length(allSpecies),
    node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
    node_end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    node_run_time = as.numeric(difftime(Sys.time(), node_start_time, units = "hours")),
    stringsAsFactors = FALSE
  )

    write.csv(log_entry, paste0(myspecies, "_log.csv"))
}

# Generate the job name with the current date
jobname <- paste0('dylcar_explore_occ_run_OCCTI_ANTS', "_", format(Sys.Date(), "%d_%m_%Y"))
dir.create(paste0("_rslurm_", jobname))

# Slurm job submission
sjob <- slurm_apply(
  f = occti_run,
  params = data.frame(species_i = 1:length(allSpecies)),
  jobname = jobname,
  nodes = length(allSpecies),
  cpus_per_node = 1,
  submit = TRUE,
  global_objects = c("allSpecies", "data", "data_range", "add_dates", "fit_occ", "fit_trend", "expit", "pcfunc", "pcfunc2", "sigfunc"),
  slurm_options = list(time = "01:00:00", mem = 20000, error = "%a.err",
  account = "ceh_generic", partition = "standard", qos = "short") ### HERE
)

#SBATCH --account=mygws
#SBATCH --partition=debug
