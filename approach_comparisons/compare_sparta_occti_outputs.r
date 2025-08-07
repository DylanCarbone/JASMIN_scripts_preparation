library(dplyr)
library(sparta)
library(ggplot2)
library(patchwork)
library(mgcv)
library(ggplot2)
library(gtools)

ants_sparta_paths_RW = list.files("past_jasmin_datalabs_runs/datalabs_ants_sparta_RW", pattern = "*.rdata", full.names = TRUE)
ants_sparta_paths = list.files("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_SPARTA_ANTS_NO_RW_21_05_2025", pattern = "*.rdata", full.names = TRUE)
ants_occti_paths = list.files("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_OCCTI_ANTS_17_04_2025/results", pattern = "*nstart_5_occupancy_output.rds", full.names = TRUE)

sparta_ants_outputs_RW = list()

# Create a new environment
tmp_env <- new.env()

for (path_i in 1:length(ants_sparta_paths_RW)){

  # Load the .RData file into that environment
  load(ants_sparta_paths_RW[path_i], envir = tmp_env)

  # Or convert everything to a list
  env_contents <- as.list(tmp_env)$out

  sparta_ants_outputs_RW[[path_i]] = env_contents

  names(sparta_ants_outputs_RW)[path_i] <- env_contents$SPP_NAME

}

sparta_ants_outputs = list()

# Create a new environment
tmp_env <- new.env()

for (path_i in 1:length(ants_sparta_paths)){

# Load the .RData file into that environment
load(ants_sparta_paths[path_i], envir = tmp_env)

# Or convert everything to a list
env_contents <- as.list(tmp_env)$out

sparta_ants_outputs[[path_i]] = env_contents

names(sparta_ants_outputs)[path_i] <- env_contents$SPP_NAME

}

occti_ants_outputs = list()

for (path_i in 1:length(ants_occti_paths)){

sp_output = readRDS(ants_occti_paths[path_i])

occti_ants_outputs[[path_i]] = sp_output

names(occti_ants_outputs)[path_i] <- sp_output$Species

}

species_intersect = intersect(intersect(names(occti_ants_outputs), names(sparta_ants_outputs_RW)), names(sparta_ants_outputs))

# filter to only include species that exist in both datasets
occti_ants_outputs <- occti_ants_outputs[names(occti_ants_outputs) %in% species_intersect]
sparta_ants_outputs_RW <- sparta_ants_outputs_RW[names(sparta_ants_outputs_RW) %in% species_intersect]
sparta_ants_outputs <- sparta_ants_outputs[names(sparta_ants_outputs) %in% species_intersect]

dir.create("occupancy_comparison_plots")

bound_occ = function(num){

  num = ifelse(num >= 1, 1-1e-6, ifelse(num <= 0, 1e-6, num))

  return(num)
}

# Loop through each species
for (species in species_intersect) {
  
  # Sparta plot
  sparta_RW_plot <- plot(sparta_ants_outputs_RW[[species]]) +
    ggtitle(paste(species, "- sparta with random walk")) +
    theme(plot.title = element_text(hjust = 0.5))

  sparta_plot <- plot(sparta_ants_outputs[[species]]) +
    ggtitle(paste(species, "- sparta without random walk")) +
    theme(plot.title = element_text(hjust = 0.5))

  # occti plot
  occ_data <- occti_ants_outputs[[species]]$Index

  # Step 1: Simulate 1000 values for each year between psiA_L and psiA_U
  n_iter <- 1000

  psiA_draws <- do.call(rbind, lapply(1:nrow(occ_data), function(i) {
    year <- occ_data$Year[i]
    mean_val <- occ_data$psiA[i]
    sd_val <- (occ_data$psiA_U[i] - occ_data$psiA_L[i]) / (2 * 1.96)

    data.frame(
      Year = year,
      Iteration = 1:n_iter,
      Simulated_psiA = rnorm(n_iter, mean = mean_val, sd = sd_val)
    )
  })) %>%
  filter(Simulated_psiA < 1, Simulated_psiA > 0)

  years_range = range(occ_data$Year)
  years = years_range[1]:years_range[2]

  loess_predictions <- psiA_draws %>%
  group_by(Iteration) %>%
  do({
    mod <- loess(Simulated_psiA ~ Year, data = ., span = 0.75)
    data.frame(Year = years, Pred = predict(mod, newdata = data.frame(Year = years)))
  }) %>%
  ungroup() %>%
  mutate(Pred = bound_occ(Pred))

  loess_summary <- loess_predictions %>%
  group_by(Year) %>%
  summarise(
    psiA_loess_mean = mean(Pred, na.rm = TRUE),
    psiA_loess_lower = quantile(Pred, 0.025, na.rm = TRUE),
    psiA_loess_upper = quantile(Pred, 0.975, na.rm = TRUE)
  )

  # Plot
  occti_plot <- ggplot(occ_data, aes(x = Year)) + 
    geom_line(aes(y = psiA), colour = "blue", alpha = 0.25) +
    geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), fill = "blue", alpha = 0.15) +
    geom_line(data = loess_summary, aes(y = psiA_loess_mean), colour = "darkred", size = 1.2) +
    geom_ribbon(data = loess_summary, aes(ymin = psiA_loess_lower, ymax = psiA_loess_upper), fill = "red", alpha = 0.15) +
    labs(x = "Year", y = "Occupancy Index") +
    theme_minimal() +
    ggtitle(paste(species, "- occti with span = 0.75")) +
    ylim(0,1)

# Combine and print
combined_plot <- sparta_RW_plot + occti_plot

# ---- Save to file ----
file_name <- paste0("occupancy_comparison_plots", "/", species, "_comparison.png")
ggsave(filename = file_name, plot = combined_plot, width = 18, height = 6, dpi = 300)
}
