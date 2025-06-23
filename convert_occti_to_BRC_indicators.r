library(dplyr)
library(devtools)
library(BRCindicators)

load_all("../occLite")

# Load all OCCTI model outputs
occti_ants_outputs <- load_occti_outputs("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_OCCTI_ANTS_01_05_2025/results")

# Generate plots and return smoothed LOESS summaries
loess_results <- smooth_occti_outputs(
  occti_outputs = occti_ants_outputs,
  save_dir = "occti_plots",
  save_plots = TRUE,
  span = 0.75,
  plot_width = 18,
  plot_height = 6,
  n_iter = 1000
)

# Note that values have been logit transformed.
bma_df = convert_for_bma(loess_results)

bma_indicator = bma(bma_df, parallel = TRUE, m.scale = "logit")

plot_indicator(indicator = bma_indicator[,'Index.M'],
               CIs = bma_indicator[,c(3,4)])
