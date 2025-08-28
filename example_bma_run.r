library(dplyr)
library(devtools)
library(BRCindicators)
library(tidyr)
library(patchwork)
library(stringr)
library(tidyr)

# use local package helpers (e.g., posterior_samp_occti, load_occti_outputs)
devtools::install_github("DylanCarbone/occLite")
library(occLite)

# Trend start year
trend_start_year = 1970
trend_end_year = 2024

# I am going to load the posterior samples - They will need to be stored locally as loading them from the shared folder can take a very long time
# NB: The species outputs are all lowercase and this may affect merging with species lists. Please let me know if I should change this.
sparta_rv_beetles_all_sp_posterior_samples <- applySamp_single(
  modPath = "RoveBeetles_Sp/_rslurm_RoveBeetles",
  region = "UK",
  nSamps = 999,
  minObs = NULL,
  speciesToKeep = NA,
  t0 = trend_start_year,
  tn = trend_end_year
)

# Somewhere here we need to filter based on a species list...

# This currently does not work as species do not seem to have rule of thumb
sparta_rv_beetles_all_sp_posterior_samples <- filter_rot(sparta_rv_beetles_all_sp_posterior_samples)

# Produce an indicator
bma_sparta_rv_beetles_all_sp = calculate_indicator(sparta_rv_beetles_all_sp_posterior_samples$posterior_samples, sparta_output = TRUE, method = "bma", bma_ind = "prime", save_outputs = TRUE, outputs_path = "bma_runs/bma_sparta_rv_beetles_all_sp.rds")

bma_sparta_rv_beetles_all_sp = readRDS("bma_runs/bma_sparta_rv_beetles_all_sp.rds")

sparta_plots_rv_beetles_all_sp <- summariseMSI(bma_sparta_rv_beetles_all_sp, label = "sparta all species",
                                      minYear = trend_start_year,
                                      maxYear = trend_end_year,
                                      writePlot = TRUE)

plot = sparta_plots_rv_beetles_all_sp$indicator
ggsave("bma_plots/rv_beetles_all_sp.png", plot = plot, width = 16, height = 4, dpi = 300)