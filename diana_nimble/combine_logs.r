library(dplyr)

log_files = list.files("diana_nimble/past_JASMIN_runs", pattern = "_log.csv", recursive = TRUE, full.name = TRUE)

log_list = list()

for (i in seq_along(log_files)){

    file = log_files[i]

  log_list[[i]] = read.csv(file)
}

# Safe rbind function
bind_same_cols <- function(df_list) {
  # Find common column names across all data frames
  common_cols <- Reduce(intersect, lapply(df_list, colnames))
  
  # Subset each data frame to only the common columns
  df_list_trimmed <- lapply(df_list, function(df) df[common_cols])
  
  # Bind the trimmed data frames together
  combined_df <- do.call(rbind, df_list_trimmed)
  
  return(combined_df)
}

# Correct way to use safe_rbind()
df <- bind_same_cols(log_list)

df = df %>% select(-X) %>%
mutate(NIMBLE_initialisation_run_time = ifelse(taxa_group %in% c("Ants", "Bees"), (NIMBLE_initialisation_run_time)/60, NIMBLE_initialisation_run_time))

write.csv(df, "ant_and_bee_logs.csv")

################################

library(ggplot2)
library(dplyr)
library(tidyr)

# Compute the correct "Other_run_time"
df <- df %>%
  mutate(Other_run_time = ifelse(!is.na(NIMBLE_initialisation_run_time), node_run_time - NIMBLE_initialisation_run_time, node_run_time)
  )

# Reshape data: Now it contains "NIMBLE_initialisation_run_time" and "Other_run_time"
df_long <- df %>%
  pivot_longer(cols = c(NIMBLE_initialisation_run_time, Other_run_time),
               names_to = "run_component", values_to = "run_time")

# Create the box plot
ggplot(df_long, aes(x = taxa_group, y = run_time, fill = run_component)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Comparison of Run Times across Taxa Groups",
       x = "Taxa Group",
       y = "Run Time (seconds)",
       fill = "Run Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Compute mean percentage of initialisation time per taxa_group
percentage_runtimes <- df %>%
  group_by(taxa_group) %>%
  summarise(percentage_run_time = 100 * mean(NIMBLE_initialisation_run_time / node_run_time),
            position = mean(node_run_time) + 1) %>%
  mutate(percentage_run_time = ifelse(!is.na(percentage_run_time), percentage_run_time, 0))

# Create the stacked bar plot with correctly computed components
ggplot(df_long %>% 
         group_by(taxa_group, run_component) %>% 
         summarise(run_time = mean(run_time), .groups = "drop"), 
       aes(x = taxa_group, y = run_time, fill = run_component)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = percentage_runtimes, 
            mapping = aes(x = taxa_group, y = position, 
                          label = paste0(round(percentage_run_time, 1), "%")), 
            inherit.aes = FALSE, vjust = -0.5, size = 4) +  # Adjust vjust for correct label placement
  theme_minimal() +
  labs(x = "taxa group",
       y = "run time hours",
       fill = "run component") +
  scale_fill_manual(values = c("NIMBLE_initialisation_run_time" = "red", 
                               "Other_run_time" = "blue")) +  # Correct labels for the components
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
