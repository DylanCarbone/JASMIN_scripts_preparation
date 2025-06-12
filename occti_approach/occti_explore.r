# setup the session and R ...
module load jasr

R

if(!require("BRCmap")){
remotes::install_github("colinharrower/BRCmap")
}

library(parallel)
library(rslurm)
library(tidyverse)
library(dplyr)

library(data.table)
library(unmarked)
library(ggplot2)
library(igr)
library(BRCmap)

source("occti_approach/prepare_data.r")

for (file in list.files("occti_approach/occti_functions", full.names = TRUE)){
source(file)
}

# Model parameters
min.Recs <- 50 # 50 # number of records per species for inclusion
nyr <- 2 # minimum number of years sampled
nstart_vector <- 5 # Number of starting values

taxa_group = "Butterflies"

data = readRDS("formatted_butterfly_data.rds")
monad_country_df = read.csv("UK_grid_refs_monads.csv")
colnames(monad_country_df) = tolower(colnames(monad_country_df))

# prepare_data is a wrapper function that removes invalid grid references, removes grid references above a 1 km monad resolution, coverts references below a 1km resolution to 1km, identifies country of the monad grid reference,
data = prepare_data(data = data)

# Summarise species
speciesSummary <- data %>%
  group_by(species) %>%
  summarise(nuSiteID = length(unique(gridref)),
            nuRecs = length(species)) %>%
  arrange(desc(nuRecs))

# Define species list
allSpecies <- sort(speciesSummary$species[speciesSummary$nuRecs > min.Recs])

# Subset by the number of visits to each site, removing sites that had data recorded for only 1 year
sites_to_include = data %>% distinct(gridref, Year) %>% 
  group_by(gridref) %>%
  summarise(n_years_sampled = n()) %>%
  filter(n_years_sampled >= nyr) %>%
  pull(gridref)

length(sites_to_include)
length(unique(data$gridref))

# Filter out sites that have not been sampled over enough years
data = data %>% filter(gridref %in% sites_to_include)

regions = c(unique(data$region), "uk")

occti_run = function(species_i, region_i){
  
  node_start_time = Sys.time()

  myspecies <- allSpecies[species_i]
  print(myspecies)

  myregion = regions[region_i]
  print(myregion)

  if (myregion != "uk") {
    data_region = data %>% filter(region == myregion)
  } else {
    data_region = data
  }

  year_range = range(data_region$Year)

  data_region$northing_scaled <- as.numeric(scale(data_region$northing))
  data_region$easting_scaled  <- as.numeric(scale(data_region$easting))

  for (nstart_i in nstart_vector) {

    # Run occupancy model
    occupancy_result <- try(
      fit_occ_formatted(
        spp = myspecies,
        obdata = data_region,
        occformula = "~ northing_scaled + I(northing_scaled^2) + easting_scaled + I(easting_scaled^2)",
        detformula = "~ logLL + SEAS",
        covnames = c("northing_scaled", "easting_scaled"),
        minyear = year_range[1],
        maxyear = year_range[2],
        trendyears = year_range[1],
        nstart = nstart_i,
        engine = "C"
      ),
      silent = TRUE
    )

    # Check for error and extract message if needed
    if (inherits(occupancy_result, "try-error")) {
      status <- "failed"
      
      # Capture error message safely
      error_message <- as.character(occupancy_result)

      # Define path and save the error message
      if (!dir.exists("error_messages")) dir.create("error_messages")
      error_file <- paste0("error_messages/", myspecies, "_", myregion, "_nstart_", nstart_i, "_occti_error.txt")
      writeLines(error_message, con = error_file)
      
    } else {
      status <- "success"
    }

    # If model ran successfully, generate plot and save
    if (status == "success") {
      plot <- ggplot(occupancy_result$Index, aes(x = Year, y = psiA)) +
        geom_line(size = 1, color = "blue") +
        geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), alpha = 0.2) +
        labs(x = "Year", y = "Occupancy Index") +
        theme_minimal()

      if (!dir.exists("plots")) dir.create("plots")
      ggsave(paste0("plots/", myspecies, "_", myregion, "_nstart_", nstart_i, ".png"), plot = plot)

      if (!dir.exists("results")) dir.create("results")
      saveRDS(occupancy_result, paste0("results/", myspecies, "_", myregion, "_nstart_", nstart_i, "_occupancy_output.rds"))
    }

    # log run attributes
    log_entry <- data.frame(
      taxa_group = taxa_group,
      species_name = myspecies,
      region = myregion,
      JASMIN = TRUE,
      queue = "long-serial",
      n_nodes_requested = (length(allSpecies) * length(regions)),
      node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
      node_end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      node_run_time = as.numeric(difftime(Sys.time(), node_start_time, units = "hours")),
      n_start_val = nstart_i,
      status = status,
      stringsAsFactors = FALSE
    )

    if (!dir.exists("logs")) dir.create("logs")
    write.csv(log_entry, paste0("logs/", myspecies, "_", myregion, "_nstart_", nstart_i, "_log.csv"), row.names = FALSE)

  }
}

# Generate the job name with the current date
jobname <- paste0('explore_occ_run_OCCTI_', toupper(taxa_group), "_", format(Sys.Date(), "%d_%m_%Y"))
dir.create(paste0("_rslurm_", jobname))

# Create a dataframe with every combination of region and species for the parameters
params_df <- expand.grid(allSpecies, regions)

# Rename columns (optional)
names(params_df) <- c("species_i", "region_i")

# Slurm job submission
sjob <- slurm_apply(
  f = occti_run,
  params = params_df,
  jobname = jobname,
  nodes = nrow(params_df),
  cpus_per_node = 1,
  submit = TRUE,
  global_objects = c("allSpecies", "regions", "data", "fit_occ_formatted", "fit_trend", "expit", "pcfunc", "pcfunc2", "sigfunc", "taxa_group", "nstart_vector"),
  slurm_options = list(time = "24:00:00", mem = 30000, error = "%a.err",
  account = "ceh_generic", partition = "standard", qos = "long") ### HERE
)

#SBATCH --account=mygws
#SBATCH --partition=debug
