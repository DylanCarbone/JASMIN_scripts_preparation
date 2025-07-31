library(dplyr)
library(devtools)
library(BRCindicators)
library(tidyr)

load_all("../occLite")
# load_all("../wrappeR")

# Load all occti model outputs
occti_ants_outputs <- load_occti_outputs("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_OCCTI_ANTS_17_04_2025/results")
# occti_butterflies_outputs <- load_occti_outputs("past_jasmin_datalabs_runs/_rslurm_explore_occ_run_OCCTI_BUTTERFLIES_12_06_2025/results", region_filter = "uk")

# Generate plots and return smoothed LOESS summaries
loess_results_ants <- smooth_occti_outputs(
  occti_outputs = occti_ants_outputs,
  save_dir = "occti_plots_ants",
  save_plots = TRUE,
  span = 0.75,
  plot_width = 18,
  plot_height = 6,
  n_iter = 1000
)

# loess_results_butterflies <- smooth_occti_outputs(
#   occti_outputs = occti_ants_outputs,
#   save_dir = "occti_plots_butterflies",
#   save_plots = TRUE,
#   span = 0.75,
#   plot_width = 18,
#   plot_height = 6,
#   n_iter = 1000
# )

# Compare standard deviations of posterior samples for ants

# Attempt to load sparta outputs
model_outputs_path = list.files("past_jasmin_datalabs_runs/datalabs_ants_sparta_RW", pattern = "*.rdata", full.names = TRUE, recursive = TRUE)
species_names = gsub(".rdata","", basename(model_outputs_path))
species_names = species_names[species_names %in% names(occti_ants_outputs)]

if (!file.exists("sparta_posterior_samples")){dir.create("sparta_posterior_samples")}

result <- applySamp_single(
  modPath = "past_jasmin_datalabs_runs/datalabs_ants_sparta_RW",        # base folder
  metaPath = "not_applicable",                  # metaPath required but unused unless using "pollinators" or "priority"
  group = "Ants",                               
  region = "UK",
  indicator = "all",                            # since you're not filtering for "priority" or "pollinators"
  nSamps = 999,
  minObs = NULL,
  scaleObs = "global",
  sample = TRUE,
  write = TRUE,
  outPath = "sparta_posterior_samples",        # matches your previous output path
  speciesToKeep = NA,# convert to comma-separated string
  drop = FALSE,                                 # not used in your example, but safe to set to FALSE
  t0 = 2000,
  tn = 2020,
  clipBy = "group",
  parallel = FALSE
)

taxa_file <- list.files("sparta_posterior_samples", pattern = "all", full.names = TRUE)

# Load the file
df <- loadRfile(taxa_file) %>%
  .[[1]] %>%
  dplyr::mutate(
    tax_grp = gsub("\\_.*", "", basename(taxa_file)),
    species = tolower(gsub(" ", "_", species))
  )

metadata <- loadRfile(taxa_file) %>%
  .[[2]] %>%
  # Standardise GB/UK column names
  setNames(gsub("_GB", "_UK", names(.))) %>%
  dplyr::mutate(
    tax_grp = gsub("\\_.*", "", basename(taxa_file))
  ) %>%
  dplyr::mutate(
    species = tolower(gsub(" ", "_", .[[1]])),
    has_rot_meta = !is.na(rot_EqualWt_r_UK)
  )

# Comparison of posterior samples
str(result$samp_post)
str(loess_results_ants$loess_predictions)

str(translate_sparta(df))
str(loess_results_ants$loess_predictions)

unique(df$species)
unique(loess_results_ants$loess_predictions$species)

sparta_indicator = calculate_indicator(df, sparta_output = TRUE, method = "lambda", bma_ind = "prime")
occti_indicator = calculate_indicator(loess_results_ants$loess_predictions, sparta_output = FALSE, method = "lambda", bma_ind = "prime")

str(loess_predictions)
range(loess_predictions$pred)

sparta_plot <- summariseMSI(sparta_indicator,
                                      plotType =  "indicator", 
                                      label = "test",
                                      method = "lambda",
                                      minYear = 1980,
                                      maxYear = 2022,
                                      writePlot = TRUE,
                                      plot = TRUE)

occti_plot <- summariseMSI(occti_indicator,
                                      plotType =  "indicator", 
                                      label = "test",
                                      method = "lambda",
                                      minYear = 1980,
                                      maxYear = 2022,
                                      writePlot = TRUE,
                                      plot = TRUE)
