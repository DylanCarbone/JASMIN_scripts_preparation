library(dplyr)
library(sparta)
library(ggplot2)
library(patchwork)  # For side-by-side plots
library(mgcv)
library(ggplot2)

ants_sparta_paths = list.files("occ_comparison/ants_occ_outputs", pattern = "*.rdata", full.names = TRUE)
ants_occti_paths = list.files("occ_comparison/_rslurm_dylcar_explore_occ_run_OCCTI_ANTS_01_05_2025/results", pattern = "*nstart_5_occupancy_output.rds", full.names = TRUE)

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

names(occti_ants_outputs)[path_i] <- sp_output$species

}

# filter to only include species that exist in both datasets
occti_ants_outputs <- occti_ants_outputs[names(occti_ants_outputs) %in% names(sparta_ants_outputs)]
sparta_ants_outputs <- sparta_ants_outputs[names(sparta_ants_outputs) %in% names(occti_ants_outputs)]
length(occti_ants_outputs) == length(sparta_ants_outputs)

dir.create("occupancy_comparison_plots")

# Loop through each species
for (species in names(occti_ants_outputs)) {
  
  # Sparta plot
  sparta_plot <- plot(sparta_ants_outputs[[species]]) +
    ggtitle(paste(species, "- SPARTA")) +
    theme(plot.title = element_text(hjust = 0.5))

  # occti plot
  occ_data <- occti_ants_outputs[[species]]$Index

  loess_fit <- loess(psiA ~ Year, data = occ_data, span = .5)
  occ_data$psiA_loess <- predict(loess_fit)

  ggplot(occ_data, aes(x = Year)) + 
    geom_line(aes(y = psiA), colour = "blue", alpha = 0.25) +
    geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), fill = "blue", alpha = 0.15) +
    geom_line(aes(y = psiA_loess), colour = "forestgreen", size = 1.2) +
    labs(x = "Year", y = "Occupancy Index") +
    theme_minimal()

  # Step 1: Estimate standard errors from CI
  occ_data$psiA_se <- (occ_data$psiA_U - occ_data$psiA_L) / (2 * 1.96)

  # Step 2: Compute weights as inverse variance
  occ_data$psiA_weight <- 1 / (occ_data$psiA_se^2)

  # Step 3: Fit weighted LOESS
  loess_weighted <- loess(psiA ~ Year, data = occ_data, span = 0.25, weights = psiA_weight)
  occ_data$psiA_loess_weighted <- predict(loess_weighted)

  # Plot
  occti_plot <- ggplot(occ_data, aes(x = Year)) + 
    geom_line(aes(y = psiA), colour = "blue", alpha = 0.25) +
    geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), fill = "blue", alpha = 0.15) +
    geom_line(aes(y = psiA_loess_weighted), colour = "forestgreen", size = 1.2) +
    labs(x = "Year", y = "Occupancy Index") +
    theme_minimal()
    
    # Combine and print
    combined_plot <- sparta_plot + occti_plot

  # ---- Save to file ----
  file_name <- paste0("occupancy_comparison_plots", "/", species, "_comparison.png")
  ggsave(filename = file_name, plot = combined_plot, width = 12, height = 6, dpi = 300)
}
