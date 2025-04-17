library(dplyr)
library(sparta)
library(ggplot2)
library(patchwork)  # For side-by-side plots

ants_sparta_paths = list.files("occ_comparison/ants_occ_outputs", pattern = "*.rdata", full.names = TRUE)
ants_occti_paths = list.files("occ_comparison/_rslurm_dylcar_explore_occ_run_OCCTI_ANTS_03_04_2025/results", pattern = "*nstart_5_occupancy_output.rds", full.names = TRUE)

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

# filter to only include species that exist in both datasets
occti_ants_outputs <- occti_ants_outputs[names(occti_ants_outputs) %in% names(sparta_ants_outputs)]
sparta_ants_outputs <- sparta_ants_outputs[names(sparta_ants_outputs) %in% names(occti_ants_outputs)]
length(occti_ants_outputs) == length(sparta_ants_outputs)

# List of shared species
species_list <- names(sparta_ants_outputs)

dir.create("occupancy_comparison_plots")

# Loop through each species
for (species in names(occti_ants_outputs)) {
  
  # Sparta plot
  sparta_plot <- plot(sparta_ants_outputs[[species]]) +
    ggtitle(paste(species, "- SPARTA")) +
    theme(plot.title = element_text(hjust = 0.5))

  # occti plot
  occti_data <- occti_ants_outputs[[species]]$Index

  occti_plot <- ggplot(occti_data, aes(x = Year, y = psiA)) + 
    geom_line(size = 1, color = "blue") +
    geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), alpha = 0.2) +
    labs(x = "Year", y = "Occupancy Index") +
    ggtitle(paste(species, "- occti")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Combine and print
  combined_plot <- sparta_plot + occti_plot

  # ---- Save to file ----
  file_name <- paste0("occupancy_comparison_plots", "/", species, "_comparison.png")
  ggsave(filename = file_name, plot = combined_plot, width = 12, height = 6, dpi = 300)
}

saveRDS(occti_ants_outputs, "occti_ants_outputs_all_species.RDS")
