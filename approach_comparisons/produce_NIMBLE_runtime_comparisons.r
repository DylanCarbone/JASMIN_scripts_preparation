library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Safe rbind function
load_and_bind_same_cols <- function(df_file_vec) {

  df_list = list()

  for (df_path_i in seq_along(df_file_vec)){

    df = read.csv(df_file_vec[df_path_i])

    df_list[[df_path_i]] = df
  }
  
  # Find common column names across all data frames
  common_cols <- Reduce(intersect, lapply(df_list, colnames))
  
  # Subset each data frame to only the common columns
  df_list_trimmed <- lapply(df_list, function(df) df[common_cols])
  
  # Bind the trimmed data frames together
  combined_df <- do.call(rbind, df_list_trimmed)
  
  return(combined_df)
}

sparta_run = load_and_bind_same_cols(list.files("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_SPARTA_ANTS_22_07_2025", pattern = "*.csv", recursive = TRUE, full.names = TRUE))

nimble_standard = load_and_bind_same_cols(list.files("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_DIANA_NIMBLE_ANTS_24_01_2025", pattern = "*.csv", recursive = TRUE, full.names = TRUE))

nimble_par_threads = load_and_bind_same_cols(list.files("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_DIANA_NIMBLE_ants_par_threads_22_07_2025", pattern = "*.csv", recursive = TRUE, full.names = TRUE))

nimble_par_species = load_and_bind_same_cols(list.files("past_jasmin_datalabs_runs/dylcar_explore_occ_run_DIANA_NIMBLE_ANTS_par_species_13_02_2025", pattern = "*.csv", recursive = TRUE, full.names = TRUE))

nimble_bees_standard = load_and_bind_same_cols(list.files("past_jasmin_datalabs_runs/_rslurm_dylcar_explore_occ_run_DIANA_NIMBLE_BEES_27_01_2025", pattern = "*.csv", recursive = TRUE, full.names = TRUE))

sparta_run_formatted = sparta_run %>%
rename(qos = queue) %>%
select(-n_start_val) %>%
mutate(job_description = "Ants with sparta", queue = "standard")

nimble_standard_formatted = nimble_standard %>%
rename(qos = queue, config_time = NIMBLE_initialisation_run_time, job_description = taxa_group) %>%
mutate(config_time = config_time/60, qos = "standard", job_description = "Ants with no parallel threading or multiple species per node (standard)") %>%
filter(node_run_time >= .5)

nimble_par_threads_formatted = nimble_par_threads %>%
rename(qos = queue, config_time = NIMBLE_initialisation_run_time, job_description = taxa_group) %>%
mutate(qos = "high", job_description = "Ants with MCMC chains handled by separate cores (parallel threads)")

nimble_par_species_formatted = nimble_par_species %>%
rename(qos = queue, config_time = NIMBLE_initialisation_run_time, job_description = taxa_group) %>%
mutate(qos = "high", job_description = "Ants with 15 species per node (parallel species)")

nimble_bees_standard_formatted = nimble_bees_standard %>%
rename(qos = queue, config_time = NIMBLE_initialisation_run_time, job_description = taxa_group) %>%
mutate(qos = "standard", job_description = "Bees with no parallel threading or multiple species per node (standard)", config_time = config_time/60)

str(sparta_run_formatted)
str(nimble_standard_formatted)
str(nimble_par_threads_formatted)
str(nimble_par_species_formatted)
str(nimble_bees_standard_formatted)

###############
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

# ---- Plot 1: Mean node_run_time with SD bars for Sparta vs Nimble Standard ----

df1 <- bind_rows(
  sparta_run_formatted %>%
    mutate(group = "Sparta (ants)"),
  nimble_standard_formatted %>%
    mutate(group = "Standard (ants)")
)

df1_summary <- df1 %>%
  group_by(group) %>%
  summarise(
    mean_runtime = mean(node_run_time, na.rm = TRUE),
    sd_runtime = sd(node_run_time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(group = factor(group, levels = c("Sparta (ants)", "Standard (ants)")))

p1 <- ggplot(df1_summary, aes(x = group, y = mean_runtime)) +
  geom_bar(stat = "identity", width = 0.5, fill = "grey80") +
  geom_errorbar(aes(ymin = mean_runtime - sd_runtime, ymax = mean_runtime + sd_runtime), width = 0.2) +
  labs(
    x = "Job description",
    y = "Node completion time (hrs)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(margin = margin(t = 10)),
    plot.title = element_blank()
  )

ggsave("approach_comparisons/plot1_mean_runtime.png", plot = p1, width = 7, height = 5, dpi = 300, bg = "white")

# ---- Plot 2: Stacked Bar of Config Time vs Run Time ----

df2 <- bind_rows(
  nimble_standard_formatted %>% mutate(group = "Standard (ants)"),
  nimble_par_threads_formatted %>% mutate(group = "Parallel threads (ants)"),
  nimble_par_species_formatted %>% mutate(group = "Parallel species (ants)"),
  nimble_bees_standard_formatted %>% mutate(group = "Standard (bees)")
)

df2_summary <- df2 %>%
  group_by(group) %>%
  summarise(
    mean_config = mean(config_time, na.rm = TRUE),
    mean_runtime = mean(node_run_time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    compute_time = mean_runtime - mean_config,
    config_percent = round((mean_config / mean_runtime) * 100, 1),
    group = factor(group, levels = c("Standard (ants)", "Parallel threads (ants)", "Parallel species (ants)", "Standard (bees)"))
  )

df2_long <- df2_summary %>%
  select(group, config_time = mean_config, compute_time) %>%
  pivot_longer(cols = c("config_time", "compute_time"), names_to = "time_type", values_to = "hours")

p2 <- ggplot(df2_long, aes(x = group, y = hours, fill = time_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(
    data = df2_summary,
    aes(x = group, y = mean_config + 0.1, label = paste0(config_percent, "%")),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c("config_time" = "steelblue", "compute_time" = "grey80"),
    labels = c(
      "config_time" = "Model configuration time",
      "compute_time" = "Model run time"
    )
  ) +
  labs(
    x = "Job description",
    y = "Node completion time (hrs)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10))
  )

ggsave("approach_comparisons/plot2_runtime_breakdown.png", plot = p2, width = 8, height = 5, dpi = 300, bg = "white")
