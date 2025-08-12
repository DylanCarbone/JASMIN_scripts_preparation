library(dplyr)
library(devtools)
library(BRCindicators)
library(tidyr)

# load_all("../occLite")
install_github("DylanCarbone/occLite")

library(occLite)

# Load all occti model outputs
occti_ants_outputs <- load_occti_outputs("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_OCCTI_ANTS_17_04_2025/results")

# Apply a species-year filter
occti_ants_outputs_filtered <- filter_occupancy(occti_ants_outputs, ci_width_threshold = 0.8, z_score_threshold = 3)

# Generate plots and return smoothed LOESS summaries
loess_results_ants <- posterior_samp_occti(
  occti_outputs = occti_ants_outputs_filtered,
  save_dir = "occti_plots_ants_with_ci_z_score_filtering",
  save_plots = TRUE,
  span = 0.75,
  plot_width = 18,
  plot_height = 6,
  n_iter = 1000
)

# Attempt to load sparta outputs
dir.create("sparta_posterior_samples")

result <- applySamp_single(
  modPath = "past_jasmin_datalabs_runs/datalabs_ants_sparta_RW",        # base folder
  metaPath = "not_applicable",
  group = "Ants",                               
  region = "UK",
  indicator = "all",
  nSamps = 999,
  minObs = NULL,
  scaleObs = "global",
  sample = TRUE,
  write = TRUE,
  outPath = "sparta_posterior_samples",        # matches your previous output path
  speciesToKeep = NA,
  drop = FALSE,                                 # not used in your example, but safe to set to FALSE
  t0 = 2000,
  tn = 2020,
  clipBy = "group",
  parallel = FALSE
)

taxa_file <- list.files("sparta_posterior_samples", pattern = "all", full.names = TRUE)

# Load the file
df_sparta <- loadRData(taxa_file) %>%
  .[[1]] %>%
  dplyr::mutate(
    tax_grp = gsub("\\_.*", "", basename(taxa_file)),
    species = tolower(gsub(" ", "_", species))
  )

sparta_indicator = calculate_indicator(df_sparta, sparta_output = TRUE, method = "bma", bma_ind = "prime")
occti_indicator = calculate_indicator(loess_results_ants, sparta_output = FALSE, method = "bma", bma_ind = "prime")

sparta_plot <- summariseMSI(sparta_indicator,
                                      plotType =  "indicator", 
                                      label = "sparta test",
                                      method = "bma",
                                      minYear = 1980,
                                      maxYear = 2022,
                                      writePlot = TRUE,
                                      plot = TRUE)

dev.new()

occti_plot <- summariseMSI(occti_indicator,
                                      plotType =  "indicator", 
                                      label = "occti test",
                                      method = "bma",
                                      minYear = 1980,
                                      maxYear = 2022,
                                      writePlot = TRUE,
                                      plot = TRUE)