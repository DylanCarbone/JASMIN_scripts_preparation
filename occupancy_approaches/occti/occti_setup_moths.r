# setup the session and R ...
module load jasr

R

# load occLite package, or use load_all() if you have a development version
# NB: it is recommended that you clean install occLite if you have not used it for some time, as it is currently in development. The installation of dependancies can take some time. However, next time you install the package, it should be shorter as it should detect that the dependancies already are installed
remotes::install_github("DylanCarbone/occLite")

library(occLite)
library(rslurm)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(lubridate)

###
setwd("/home/users/dylcar/JASMIN_scripts_preparation/moths_sp/")
data = read.csv("NMRS_data_1970to2020_180825_cleaned_monads.csv") %>%
mutate(species = gsub("/", "_", species), year = year(ymd(date)))

###

# Set the working directory if you have not already
setwd("/gws/nopw/j04/ceh_generic/dylcar/indicators_2025/moths_sp")

# Model parameters
min.Recs <- 50 # 50 # number of records per species for inclusion
nyr <- 2 # minimum number of years sampled
n_start <- 5 # Number of starting valuesCheesit9

# specify group
taxa_group = "Moths"

# Formatting required for NEW butterfly dataset only
data = read.csv("NMRS_data_1970to2020_180825.csv")

data = data %>% rename(species = taxon, date = record_date, gridref = km_sq) %>% 
select(species, date, gridref) %>%
mutate(date = dmy(date), species = gsub("/", "_", species), year = year(date))

# Prepare_data is a wrapper function that removes invalid grid references, removes grid references above a 1 km monad resolution, coverts references below a 1km resolution to 1km, identifies country of the monad grid reference,
data = prep_occ_data(data = data, subset = TRUE, min.Recs = min.Recs, nyr = nyr)

#save the data so we have a record of it post-formatting
write.csv(data, "NMRS_data_1970to2020_180825_cleaned_monads.csv")

# Obtain all species after subsetting
allSpecies = unique(data$species)

# Obtain all regions after subsetting
regions = c(unique(data$region), "gb", "uk")

occti_run = function(species, region_name){

  print(species)
  print(region_name)

  node_start_time = Sys.time()

  # Filtering for regions
  if (region_name == "uk"){
    data_region = data
  } else if (region_name == "gb"){
    data_region = data %>% filter(region != "ni")
  } else {
    data_region = data %>% filter(region == region_name)
  }
  print(paste("number of rows in region subset:", nrow(data_region)))

  # Date range
  year_range = unique(data_region$Year)

  # Scaled coords
  data_region$northing_scaled <- as.numeric(scale(data_region$northing))
  data_region$easting_scaled  <- as.numeric(scale(data_region$easting))

  # --- helper: run once and classify success/failure ---
  run_once <- function(attempt) {

    res <- try(
      fit_occ_formatted(
        spp = species,
        obdata = data_region,
        occformula = "~ northing_scaled + I(northing_scaled^2) + easting_scaled + I(easting_scaled^2)",
        detformula = "~ logLL + SEAS",
        covnames = c("northing_scaled", "easting_scaled"),
        minyear = year_range[1],
        maxyear = year_range[2],
        trendyears = year_range[1],
        nstart = n_start,
        engine = "C"
      ),
      silent = TRUE
    )

    # classify
    if (inherits(res, "try-error")) {
      list(completed = FALSE, result = NULL, msg = as.character(res))
    } else if (is.null(res$Index) || nrow(res$Index) == 0) {
      list(completed = FALSE, result = NULL, msg = paste0("empty Index on attempt ", attempt))
    } else {
      list(completed = TRUE,  result = res,  msg = NULL)
    }
  }

  # --- retry up to 2 times (3 total attempts) ---
  max_attempts <- 3
  status <- "failed"
  last_msg <- NULL
  occupancy_result <- NULL

  for (attempt in 1:max_attempts) {
    out <- run_once(attempt)
    if (out$completed) {
      status <- "success"
      occupancy_result <- out$result
      break
    } else {
      last_msg <- out$msg
      message(out$msg)
    }
  }

  # If success: plot + save results
  if (status == "success") {

    if (!dir.exists("plots")) dir.create("plots")
    if (!dir.exists("results")) dir.create("results")

    p <- ggplot2::ggplot(occupancy_result$Index, ggplot2::aes(x = Year, y = psiA)) +
      ggplot2::geom_line(size = 1, color = "blue") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = psiA_L, ymax = psiA_U), alpha = 0.2) +
      ggplot2::labs(x = "Year", y = "Occupancy Index") +
      ggplot2::theme_minimal()

    ggplot2::ggsave(paste0("plots/", species, "_", region_name, ".png"), plot = p)

    saveRDS(occupancy_result, paste0("results/", species, "_", region_name, "_occupancy_output.rds"))

  } else {
    # All attempts failed: save final error message
    if (!dir.exists("error_messages")) dir.create("error_messages")
    error_file <- paste0("error_messages/", species, "_", region_name, "_occti_error.txt")
    msg <- if (is.null(last_msg)) "Unknown error (no message captured)" else last_msg
    writeLines(msg, con = error_file)
  }

  # Log
  if (!dir.exists("logs")) dir.create("logs")
  log_entry <- data.frame(
    taxa_group = taxa_group,
    species_name = species,
    region = region_name,
    JASMIN = TRUE,
    queue = "long-serial",
    n_nodes_requested = (length(allSpecies) * length(regions)),
    node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
    node_end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    node_run_time = as.numeric(difftime(Sys.time(), node_start_time, units = "hours")),
    n_start_val = n_start,
    status = status,
    stringsAsFactors = FALSE
  )
  write.csv(log_entry, paste0("logs/", species, "_", region_name, "_log.csv"), row.names = FALSE)
}

# Generate the job name with the current date
jobname <- "explore_occ_run_OCCTI_MOTHS_22_08_2025"# paste0('explore_occ_run_OCCTI_', toupper(taxa_group), "_", format(Sys.Date(), "%d_%m_%Y"))

# Create a dataframe with every combination of region and species for the parameters
params_df <- expand.grid(allSpecies, unique(data$region))
gb_uk_df <- expand.grid(allSpecies, c("gb", "uk"))

# Rename columns (optional)
names(params_df) <- c("species", "region_name")
names(gb_uk_df) <- c("species", "region_name")

# Keep only species-region combinations in params_df that are present in data
params_df <- params_df %>%
  semi_join(distinct(data, species, region), by = c("species" = "species", "region_name" = "region")) %>%
  rbind(gb_uk_df)

# Exclude species that have already ran
sp_region_skip = list.files(file.path(paste0("_rslurm_", jobname), "results"), pattern = "*.rds")
params_df = params_df %>% filter(!paste0(species, "_", region_name, "_occupancy_output.rds") %in% sp_region_skip)
print(params_df)

# NB: occLite needs to be installed locally with install_github to allow for the nodes to access functions. Nodes will access the functions as they are in the installed version, not in the state they are in after you called load_all()

# Slurm job submission
sjob <- slurm_apply(
  f = occti_run,
  params = params_df,
  jobname = jobname,
  nodes = nrow(params_df),
  cpus_per_node = 1,
  submit = TRUE,
  global_objects = c("allSpecies", "regions", "data", "taxa_group", "n_start"),
  slurm_options = list(time = "02-00:00:00", mem = 30000, error = "%a.err",
  account = "ceh_generic", partition = "standard", qos = "long")
)
