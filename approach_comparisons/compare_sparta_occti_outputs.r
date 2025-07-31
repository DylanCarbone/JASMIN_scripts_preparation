library(dplyr)
library(sparta)
library(ggplot2)
library(patchwork)
library(mgcv)
library(ggplot2)
library(gtools)

ants_sparta_paths_RW = list.files("past_jasmin_datalabs_runs/datalabs_ants_sparta_RW", pattern = "*.rdata", full.names = TRUE)
ants_sparta_paths = list.files("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_SPARTA_ANTS_NO_RW_21_05_2025", pattern = "*.rdata", full.names = TRUE)
ants_occti_paths = list.files("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_OCCTI_ANTS_01_05_2025/results", pattern = "*nstart_5_occupancy_output.rds", full.names = TRUE)

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

  num = ifelse(num == 1, 1-1e-6, ifelse(num == 0, 1e-6, num))
  
  if(any(num > 1 | num < 0)){
    stop("Values are greater than 1 or less than 0")
  }

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
  }))

  loess_predictions <- psiA_draws %>%
  group_by(Iteration) %>%
  do({
    mod <- try(loess(Simulated_psiA ~ Year, data = ., span = 0.50), silent = TRUE)
    data.frame(Year = .$Year, Pred = predict(mod, newdata = data.frame(Year = .$Year)))
  }) %>%
  ungroup()

  loess_summary <- loess_predictions %>%
  group_by(Year) %>%
  summarise(
    psiA_loess_mean = mean(Pred, na.rm = TRUE),
    psiA_loess_lower = quantile(Pred, 0.025, na.rm = TRUE),
    psiA_loess_upper = quantile(Pred, 0.975, na.rm = TRUE)
  )

  # Plot
  occti_plot_1 <- ggplot(occ_data, aes(x = Year)) + 
    geom_line(aes(y = psiA), colour = "blue", alpha = 0.25) +
    geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), fill = "blue", alpha = 0.15) +
    geom_line(data = loess_summary, aes(y = psiA_loess_mean), colour = "darkred", size = 1.2) +
    geom_ribbon(data = loess_summary, aes(ymin = psiA_loess_lower, ymax = psiA_loess_upper), fill = "red", alpha = 0.15) +
    labs(x = "Year", y = "Occupancy Index") +
    theme_minimal() +
    ggtitle(paste(species, "- occti with span = 0.50")) +
    ylim(0,1)

  loess_predictions <- psiA_draws %>%
  group_by(Iteration) %>%
  do({
    mod <- try(loess(Simulated_psiA ~ Year, data = ., span = 0.75), silent = TRUE)
    data.frame(Year = .$Year, Pred = predict(mod, newdata = data.frame(Year = .$Year)))
  }) %>%
  ungroup()

  loess_summary <- loess_predictions %>%
  group_by(Year) %>%
  summarise(
    psiA_loess_mean = mean(Pred, na.rm = TRUE),
    psiA_loess_lower = quantile(Pred, 0.025, na.rm = TRUE),
    psiA_loess_upper = quantile(Pred, 0.975, na.rm = TRUE)
  )

  # Plot
  occti_plot_2 <- ggplot(occ_data, aes(x = Year)) + 
    geom_line(aes(y = psiA), colour = "blue", alpha = 0.25) +
    geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), fill = "blue", alpha = 0.15) +
    geom_line(data = loess_summary, aes(y = psiA_loess_mean), colour = "darkred", size = 1.2) +
    geom_ribbon(data = loess_summary, aes(ymin = psiA_loess_lower, ymax = psiA_loess_upper), fill = "red", alpha = 0.15) +
    labs(x = "Year", y = "Occupancy Index") +
    theme_minimal() +
    ggtitle(paste(species, "- occti with span = 0.75")) +
    ylim(0,1)

# # Simulate from log-transformed occupancy index
# psiA_draws <- do.call(rbind, lapply(1:nrow(occ_data), function(i) {
#   year <- occ_data$Year[i]
#   mean_val <- occ_data$psiA[i]
#   sd_val <- (occ_data$psiA_U[i] - occ_data$psiA_L[i]) / (2 * 1.96)

#   sim_vals <- rnorm(n_iter, mean = mean_val, sd = sd_val)
#   sim_vals <- ifelse(sim_vals > 0, sim_vals, 1e-6)
#   sim_vals <- ifelse(sim_vals > 1, 1, sim_vals)
#   sim_vals <- log(sim_vals)

#   data.frame(
#     Year = year,
#     Iteration = 1:n_iter,
#     Simulated_psiA = sim_vals
#   )
# }))

# # Apply LOESS smoothing to each iteration
# loess_predictions <- psiA_draws %>%
#   group_by(Iteration) %>%
#   do({
#     mod <- try(loess(Simulated_psiA ~ Year, data = ., span = 0.75), silent = TRUE)
#     if (inherits(mod, "try-error")) {
#       return(data.frame(Year = .$Year, Pred = NA))
#     }
#     data.frame(Year = .$Year, Pred = predict(mod, newdata = data.frame(Year = .$Year)))
#   }) %>%
#   ungroup()

# # Summarise across iterations using quantiles
# loess_summary <- loess_predictions %>%
#   group_by(Year) %>%
#   summarise(
#     psiA_loess_mean = mean(Pred, na.rm = TRUE),
#     psiA_loess_lower = quantile(Pred, 0.025, na.rm = TRUE),
#     psiA_loess_upper = quantile(Pred, 0.975, na.rm = TRUE)
#   )

# occti_plot_log <- ggplot() +
#   # Raw log-transformed occupancy index
#   geom_line(data = occ_data, aes(x = Year, y = log(psiA)), colour = "blue", alpha = 0.25) +
#   geom_ribbon(data = occ_data, aes(x = Year, ymin = log(psiA_L), ymax = log(psiA_U)), fill = "blue", alpha = 0.15) +

#   # LOESS smoothed curve and 95% interval from draws
#   geom_line(data = loess_summary, aes(x = Year, y = psiA_loess_mean), colour = "darkred", size = 1.2) +
#   geom_ribbon(data = loess_summary, aes(x = Year, ymin = psiA_loess_lower, ymax = psiA_loess_upper), fill = "red", alpha = 0.15) +

#   labs(x = "Year", y = "log(Occupancy Index)") +
#   theme_minimal() +
#   ggtitle(paste(species, "- log-transformed occupancy with 95% percentile"))

# # Guard functions
# safe_exp <- function(x) pmin(pmax(exp(x), 0), 1)  # Avoid values > 1 if smoothing overshoots

# # Summarise and back-transform
# loess_summary <- loess_predictions %>%
#   group_by(Year) %>%
#   summarise(
#     psiA_loess_mean = safe_exp(mean(Pred, na.rm = TRUE)),
#     psiA_loess_lower = safe_exp(quantile(Pred, 0.025, na.rm = TRUE)),
#     psiA_loess_upper = safe_exp(quantile(Pred, 0.975, na.rm = TRUE))
#   )

# # Plot on original [0, 1] scale (back-transformed)
# occti_plot_log_backtransformed <- ggplot() +
#   # Raw occupancy values
#   geom_line(data = occ_data, aes(x = Year, y = psiA), colour = "blue", alpha = 0.25) +
#   geom_ribbon(data = occ_data, aes(x = Year, ymin = psiA_L, ymax = psiA_U), fill = "blue", alpha = 0.15) +

#   # Smoothed mean and quantiles from log-scale back-transformed
#   geom_line(data = loess_summary, aes(x = Year, y = psiA_loess_mean), colour = "darkred", size = 1.2) +
#   geom_ribbon(data = loess_summary, aes(x = Year, ymin = psiA_loess_lower, ymax = psiA_loess_upper), fill = "red", alpha = 0.15) +

#   labs(x = "Year", y = "Occupancy Index") +
#   theme_minimal() +
#   ggtitle(paste(species, "- occupancy with 95% quantile (log back-transformed)")) +
#   ylim(0,1)

occ_data_logit = occ_data %>%
mutate(psiA = logit(bound_occ(psiA)),
psiA_L = logit(bound_occ(psiA_L)), 
psiA_U = logit(bound_occ(psiA_U)))

# Simulate from log-transformed occupancy index
psiA_draws <- do.call(rbind, lapply(1:nrow(occ_data_logit), function(i) {
  year <- occ_data_logit$Year[i]
  mean_val <- (occ_data_logit$psiA[i])
  sd_val <- (occ_data_logit$psiA_U[i] - occ_data_logit$psiA_L[i]) / (2 * 1.96)

  sim_vals <- rnorm(n_iter, mean = mean_val, sd = sd_val)

  data.frame(
    Year = year,
    Iteration = 1:n_iter,
    Simulated_psiA = sim_vals
  )
}))

# Apply LOESS smoothing to each iteration
loess_predictions <- psiA_draws %>%
  group_by(Iteration) %>%
  do({
    mod <- try(loess(Simulated_psiA ~ Year, data = ., span = 0.75), silent = TRUE)
    if (inherits(mod, "try-error")) {
      return(data.frame(Year = .$Year, Pred = NA))
    }
    data.frame(Year = .$Year, Pred = predict(mod, newdata = data.frame(Year = .$Year)))
  }) %>%
  ungroup()

# Summarise across iterations using quantiles
loess_summary <- loess_predictions %>%
  group_by(Year) %>%
  summarise(
    psiA_loess_mean = mean(Pred, na.rm = TRUE),
    psiA_loess_lower = quantile(Pred, 0.025, na.rm = TRUE),
    psiA_loess_upper = quantile(Pred, 0.975, na.rm = TRUE)
  )

occti_plot_logit <- ggplot() +

  # LOESS smoothed curve and 95% interval from draws
  geom_line(data = loess_summary, aes(x = Year, y = psiA_loess_mean), colour = "darkred", size = 1.2) +
  geom_ribbon(data = loess_summary, aes(x = Year, ymin = psiA_loess_lower, ymax = psiA_loess_upper), fill = "red", alpha = 0.15) +

  labs(x = "Year", y = "log(Occupancy Index)") +
  theme_minimal() +
  ggtitle(paste(species, "- logit-transformed occupancy with 95% percentile"))

occ_data_logit_inv = occ_data_logit %>%
mutate(psiA = inv.logit(psiA),
psiA_L =  inv.logit(psiA_L), 
psiA_U =  inv.logit(psiA_U))

# Simulate from log-transformed occupancy index
psiA_draws <- do.call(rbind, lapply(1:nrow(occ_data_logit_inv), function(i) {
  year <- occ_data_logit_inv$Year[i]
  mean_val <- (occ_data_logit_inv$psiA[i])
  sd_val <- (occ_data_logit_inv$psiA_U[i] - occ_data_logit_inv$psiA_L[i]) / (2 * 1.96)

  sim_vals <- rnorm(n_iter, mean = mean_val, sd = sd_val)

  data.frame(
    Year = year,
    Iteration = 1:n_iter,
    Simulated_psiA = sim_vals
  )
}))

# Apply LOESS smoothing to each iteration
loess_predictions <- psiA_draws %>%
  group_by(Iteration) %>%
  do({
    mod <- try(loess(Simulated_psiA ~ Year, data = ., span = 0.75), silent = TRUE)
    if (inherits(mod, "try-error")) {
      return(data.frame(Year = .$Year, Pred = NA))
    }
    data.frame(Year = .$Year, Pred = predict(mod, newdata = data.frame(Year = .$Year)))
  }) %>%
  ungroup()

# Summarise across iterations using quantiles
loess_summary <- loess_predictions %>%
  group_by(Year) %>%
  summarise(
    psiA_loess_mean = mean(Pred, na.rm = TRUE),
    psiA_loess_lower = quantile(Pred, 0.025, na.rm = TRUE),
    psiA_loess_upper = quantile(Pred, 0.975, na.rm = TRUE)
  )

occti_plot_logit_inv <- ggplot() +

  # LOESS smoothed curve and 95% interval from draws
  geom_line(data = loess_summary, aes(x = Year, y = psiA_loess_mean), colour = "darkred", size = 1.2) +
  geom_ribbon(data = loess_summary, aes(x = Year, ymin = psiA_loess_lower, ymax = psiA_loess_upper), fill = "red", alpha = 0.15) +

  labs(x = "Year", y = "log(Occupancy Index)") +
  theme_minimal() +
  ggtitle(paste(species, "- backtransformed from logit occupancy with 95% percentile")) +
  ylim(0,1)

# Combine and print
combined_plot <- sparta_RW_plot + occti_plot_2 #+ sparta_plot + occti_plot_1 + occti_plot_2 + occti_plot_logit + occti_plot_logit_inv

# ---- Save to file ----
file_name <- paste0("occupancy_comparison_plots", "/", species, "_comparison.png")
ggsave(filename = file_name, plot = combined_plot, width = 18, height = 6, dpi = 300)
}
