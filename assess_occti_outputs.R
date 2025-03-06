library(dplyr)
library(rnrfa)

library(data.table)
library(unmarked)
library(ggplot2)
library(igr)

# Save the results
results_files = list.files("past_JASMIN_runs/_rslurm_dylcar_explore_occ_run_OCCTI_BUTTERFLIES_05_03_2025", pattern = ".rds", recursive = TRUE, full.name = TRUE)

occupancy_results = readRDS(results_files[1])

# Extract trend data for a species
trend_data <- occupancy_results$Index

# Plot occupancy trend over time
ggplot(trend_data, aes(x = Year, y = psiA)) +
  geom_line(size = 1, color = "blue") +
  geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), alpha = 0.2) +
  labs(title = "Occupancy Trend for Lep_5779",
       x = "Year", y = "Occupancy Index") +
  theme_minimal()

# Compute trends for "Lep_5779"
trend_results <- fit_trend(startyear = 2015, z = occupancy_results[["tik_109"]]$Index)

# Display trend estimates
print(trend_results)
