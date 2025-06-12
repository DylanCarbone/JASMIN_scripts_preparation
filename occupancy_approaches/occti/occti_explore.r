# setup the session and R ...
module load jasr

R

if(!require("BRCmap")){
remotes::install_github("colinharrower/BRCmap")
}

# load occLite package, or use load_all() if you have a development version
# if(!require("occLite")){
# remotes::install_github("DylanCarbone/occLite")
# }
devtools::load_all("../occLite")

library(rslurm)
library(dplyr)

# Model parameters
min.Recs <- 50 # 50 # number of records per species for inclusion
nyr <- 2 # minimum number of years sampled
nstart_vector <- 5 # Number of starting values

# Butterflies group
taxa_group = "Butterflies"

# load butterfly data
data = readRDS("formatted_butterfly_data.rds")

# Prepare_data is a wrapper function that removes invalid grid references, removes grid references above a 1 km monad resolution, coverts references below a 1km resolution to 1km, identifies country of the monad grid reference,
data = prep_occ_data(data = data, subset = TRUE, min.Recs = 10, nyr = 2)

# Obtain all species after subsetting
allSpecies = unique(data$species)

# Obtain all regions after subsetting
regions = c(unique(data$region), "gb", "uk")

# Prepare function to pass to each node
occti_run = function(species, region){
  
  # Record start time for later logging
  node_start_time = Sys.time()

  # Filtering for regions
  if (region == "uk"){
    data_region = data
  } else if (region == "gb"){
    data_region = data %>% filter(region != "ir")
  } else {
    data_region = data %>% filter(region == region)
  }

  # Obtain date range
  year_range = range(data_region$Year)

  # Scale northing and easting
  data_region$northing_scaled <- as.numeric(scale(data_region$northing))
  data_region$easting_scaled  <- as.numeric(scale(data_region$easting))

  # For each n_start value specified
  for (nstart_i in nstart_vector) {

    # Run occupancy model with occti
    occupancy_result <- try(
      fit_occ_formatted(
        spp = species,
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
      error_file <- paste0("error_messages/", species, "_", region, "_nstart_", nstart_i, "_occti_error.txt")
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

      # Save plots
      if (!dir.exists("plots")) dir.create("plots")
      ggsave(paste0("plots/", species, "_", region, "_nstart_", nstart_i, ".png"), plot = plot)

      # Save Results
      if (!dir.exists("results")) dir.create("results")
      saveRDS(occupancy_result, paste0("results/", species, "_", region, "_nstart_", nstart_i, "_occupancy_output.rds"))
    }

    # log run attributes
    log_entry <- data.frame(
      taxa_group = taxa_group,
      species_name = species,
      region = region,
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

    # Save log entry
    if (!dir.exists("logs")) dir.create("logs")
    write.csv(log_entry, paste0("logs/", species, "_", region, "_nstart_", nstart_i, "_log.csv"), row.names = FALSE)

  }
}

# Generate the job name with the current date
jobname <- paste0('explore_occ_run_OCCTI_', toupper(taxa_group), "_", format(Sys.Date(), "%d_%m_%Y"))
dir.create(paste0("_rslurm_", jobname))

# Create a dataframe with every combination of region and species for the parameters
params_df <- expand.grid(allSpecies, unique(data$region))
gb_uk_df <- expand.grid(allSpecies, c("gb", "uk"))

# Rename columns (optional)
names(params_df) <- c("species", "region")
names(gb_uk_df) <- c("species", "region")

# Keep only species-region combinations in params_df that are present in data
params_df <- params_df %>%
  semi_join(distinct(data, species, region), by = c("species", "region")) %>%
  rbind(gb_uk_df)

# NB: occLite needs to be installed locally to allow for the nodes to access functions.
# Likewise, the nodes will access the functions as they are in the installed version, not in the state they are in after you called load_all()

# Slurm job submission
sjob <- slurm_apply(
  f = occti_run,
  params = params_df,
  jobname = jobname,
  nodes = nrow(params_df),
  cpus_per_node = 1,
  submit = TRUE,
  global_objects = c("allSpecies", "regions", "data", "taxa_group", "nstart_vector"),
  slurm_options = list(time = "24:00:00", mem = 30000, error = "%a.err",
  account = "ceh_generic", partition = "standard", qos = "long")
)
