library(dplyr)
library(devtools)
library(BRCindicators)
library(tidyr)

load_all("../occLite")

most_common_species = readRDS("occurence_datasets/monad_occupancy_dataset_ants.rds") %>%
count(tik) %>%
arrange(desc(n)) %>%
slice(1:10) %>%
pull(tik)

# Load all occti model outputs
occti_ants_outputs <- load_occti_outputs("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_OCCTI_ANTS_17_04_2025/results")

# Generate plots and return smoothed LOESS summaries
loess_results_ants <- posterior_samp_occti(
  occti_outputs = occti_ants_outputs,
  save_dir = "occti_plots_ants",
  save_plots = FALSE,
  span = 0.75,
  plot_width = 18,
  plot_height = 6,
  n_iter = 1000
) %>%
filter(species %in% most_common_species)

# Attempt to load sparta outputs
dir.create("sparta_posterior_samples")

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
  speciesToKeep = paste(most_common_species, collapse = ","),
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

# metadata <- loadRData(taxa_file) %>%
#   .[[2]] %>%
#   # Standardise GB/UK column names
#   setNames(gsub("_GB", "_UK", names(.))) %>%
#   dplyr::mutate(
#     tax_grp = gsub("\\_.*", "", basename(taxa_file))
#   ) %>%
#   dplyr::mutate(
#     species = tolower(gsub(" ", "_", .[[1]])),
#     has_rot_meta = !is.na(rot_EqualWt_r_UK)
#   )

# Comparison of posterior samples

str(translate_sparta(df_sparta))
str(loess_results_ants)
range(loess_results_ants$pred)
range(translate_sparta(df_sparta)$pred)

boxplot(loess_results_ants$pred)
dev.new()
boxplot(translate_sparta(df_sparta)$pred)

unique(df_sparta$species)
unique(loess_results_ants$species)

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

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)

# Step 1: Combine trend estimates from both models
sparta_trends <- sparta_indicator$final$species_assessment %>%
  rownames_to_column("species") %>%
  rename(sparta_change = percent_change_year)

occti_trends <- occti_indicator$final$species_assessment %>%
  rownames_to_column("species") %>%
  rename(occti_change = percent_change_year)

comparison_df <- sparta_trends %>%
  left_join(occti_trends, by = "species")

# Step 2: Convert to long format for plotting
plot_data <- comparison_df %>%
  pivot_longer(cols = c(sparta_change, occti_change),
               names_to = "method", values_to = "change") %>%
  mutate(method_label = ifelse(method == "sparta_change", "Sparta", "OCCTI")) %>%
  mutate(start_year = 1980,
         end_year = 2022,
         start_val = 100,
         end_val = 100 + change)

# Step 3: Plot per species
species_plots <- plot_data %>%
  select(species, method_label, start_year, end_year, start_val, end_val, change) %>%
  pivot_longer(cols = c(start_year, end_year), names_to = "point", values_to = "year") %>%
  pivot_longer(cols = c(start_val, end_val), names_to = "value_type", values_to = "occupancy") %>%
  filter((point == "start_year" & value_type == "start_val") |
         (point == "end_year" & value_type == "end_val")) %>%
  select(species, method_label, year, occupancy, change) %>%
  group_by(species) %>%
  group_split() %>%
  map(~{
    sp_data <- .x
    sp <- sp_data$species[1]

    sparta_val <- round(sp_data$change[sp_data$method_label == "Sparta"][1], 1)
    occti_val  <- round(sp_data$change[sp_data$method_label == "OCCTI"][1], 1)

    ggplot(sp_data, aes(x = year, y = occupancy, colour = method_label, group = method_label)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      theme_minimal() +
      scale_colour_manual(values = c("Sparta" = "blue", "OCCTI" = "red")) +
      ylim(0, max(sp_data$occupancy, na.rm = TRUE) * 1.1) +
      labs(
        title = paste0(sp, "\nSparta: ", sparta_val, "%/yr   |   OCCTI: ", occti_val, "%/yr"),
        x = "Year",
        y = "Relative occupancy (%)",
        colour = "Method"
      )
  })

# Step 4: View or save
# species_plots[[1]]
walk2(species_plots, comparison_df$species, ~ggsave(paste0("sparta_occti_trends_plots/plot_", .y, ".png"), .x, width = 6, height = 4))
