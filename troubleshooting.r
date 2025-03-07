library(data.table)
library(unmarked)
library(parallel)
library(lubridate)
library(ggplot2)
library(dplyr)
library(rnrfa)
library(igr)

for (file in list.files("occti_functions", full.names = TRUE)){
source(file)
}

# Read and preprocess data
data <- readRDS("monad_occupancy_dataset_ants.rds") %>% 
  rename(Species = tik, SiteID = GRIDREF, Date = lower_date) %>%
  mutate(Date = as.Date(Date),
         yday = lubridate::yday(Date),
         Year = lubridate::year(Date))

# Rename columns to match expected function inputs
setnames(data, "SiteID", "Gridref")  # Rename SiteID to Gridref


gb_ir <- split(data, igr_is_valid(data$Gridref))
names(gb_ir) = ifelse(length(gb_ir) ==2, c("gb", "ir"), "gb")

# For now, I will use only british data
data = gb_ir$gb

# Convert Date column and extract Year, Month, Week
data <- add_dates(data)

# Calculate list length (number of species observed per site-date)
data <- add_listL(data)

# Summarise species
speciesSummary <- data %>%
  group_by(Species) %>%
  summarise(nuSiteID = length(unique(Gridref)),
            nuRecs = length(Species)) %>%
  arrange(desc(nuRecs))

# Define species list
allSpecies <- sort(speciesSummary$Species[speciesSummary$nuRecs > min.Recs])

combinations <- c(
  "SA", "TA", "NA", "HA", "OA", "IA", "SB", "TB", "NB", "HB", "OB", "IB", "SC", "TC", "NC",
  "HC", "OC", "IC", "SD", "TD", "ND", "HD", "OD", "ID", "SE", "TE", "NE", "HE", "OE", "IE",
  "SF", "TF", "NF", "HF", "OF", "IF", "SG", "TG", "NG", "HG", "OG", "IG", "SH", "TH", "NH",
  "HH", "OH", "IH", "SJ", "TJ", "NJ", "HJ", "OJ", "IJ", "SK", "TK", "NK", "HK", "OK", "IK",
  "SL", "TL", "NL", "HL", "OL", "IL", "SM", "TM", "NM", "HM", "OM", "IM", "SN", "TN", "NN",
  "HN", "ON", "IN", "SO", "TO", "NO", "HO", "OO", "IO", "SP", "TP", "NP", "HP", "OP", "IP",
  "SQ", "TQ", "NQ", "HQ", "OQ", "IQ", "SR", "TR", "NR", "HR", "OR", "IR", "SS", "TS", "NS",
  "HS", "OS", "IS", "ST", "TT", "NT", "HT", "OT", "IT", "SU", "TU", "NU", "HU", "OU", "IU",
  "SV", "TV", "NV", "HV", "OV", "IV", "SW", "TW", "NW", "HW", "OW", "IW", "SX", "TX", "NX",
  "HX", "OX", "IX", "SY", "TY", "NY", "HY", "OY", "IY", "SZ", "TZ", "NZ", "HZ", "OZ", "IZ"
)

data <- data[data$Gridref %in% grep(paste0("^(", paste(combinations, collapse = "|"), ")"), data$Gridref, value = TRUE), ]

positions = rnrfa::osg_parse(grid_refs = data$Gridref, coord_system = "BNG")

data$North = positions$northing
data$East = positions$easting

data_range = range(data$Year)

# Calculate list length (number of species observed per site-date)
data <- add_listL(data)

# Check structure
str(data)

# Define species list
species_list <- unique(data$Species)

#####################

# Run multi-species occupancy model
occupancy_results_ms <- fit_occ_ms(
  splist = species_list,
  obdata = data,
  occformula = "~ North + I(North^2) + East + I(East^2)",  # Quadratic spatial covariates
  detformula = "~ logLL + SEAS",  # Detection probability model
  covnames = c("East", "North"),  # Covariates
#   parallel = TRUE,  # Enable parallel processing
#   cpus = detectCores() - 2,  # Use available CPU cores (leaving some free)
  minyear = data_range[1],
  maxyear = data_range[2],
  trendyears = data_range,
  nstart = 1,  # Number of initial values
  engine = "C"
)

# Extract trend data for a species
trend_data <- occupancy_results_ms[["tik_109"]]$Index

# Plot occupancy trend over time
ggplot(trend_data, aes(x = Year, y = psiA)) +
  geom_line(size = 1, color = "blue") +
  geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), alpha = 0.2) +
  labs(title = NA,
       x = "Year", y = "Occupancy Index") +
  theme_minimal()

# Compute trends for "Lep_5779"
trend_results <- fit_trend(startyear = 2015, z = trend_data)

# Display trend estimates
print(trend_results)


####################################

# Run the single-species occupancy model
occupancy_results_single  <- fit_occ(
spp = "tik_109",
obdata = data,
occformula = "~ North + I(North^2) + East + I(East^2)",  # Quadratic spatial covariates
detformula = "~ logLL + SEAS",  # Detection probability model
covnames = c("East", "North"),  # Covariates used
minyear = data_range[1],  # Start year
maxyear = data_range[2],  # End year
trendyears = data_range,  # Years for trend estimation
nstart = 1,  # Number of initial values to test
outputdir = "",  # Directory to save results
engine = "C"
)

failed_tally = 0

for (file in list.files("past_jasmin_runs/_rslurm_dylcar_explore_occ_run_OCCTI_BUTTERFLIES_07_03_2025", full.names = TRUE, pattern = "_occupancy_output.rds")){

  occupancy_results_single = readRDS(file)

  if(!is.null(occupancy_results_single$Index)){

  # Extract trend data for a species
  trend_data <- occupancy_results_single$Index

  # Plot occupancy trend over time
  plot = ggplot(trend_data, aes(x = Year, y = psiA)) +
    geom_line(size = 1, color = "blue") +
    geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), alpha = 0.2) +
    labs(x = "Year", y = "Occupancy Index") +
  theme_classic()  # Minimalist white background

  ggsave(paste0("plots/", occupancy_results_single$Species, ".png"), plot = plot)
}

else{ failed_tally = failed_tally + 1}

}

print(failed_tally)



# # Compute trends for "Lep_5779"
# trend_results <- fit_trend(startyear = 2015, z = trend_data)

# # Display trend estimates
# print(trend_results)



###########################


occupancy_results_single = readRDS("tik_109_occupancy_output.rds")

# Extract trend data for a species
trend_data <- occupancy_results_single$Index

# Plot occupancy trend over time
ggplot(trend_data, aes(x = Year, y = psiA)) +
  geom_line(size = 1, color = "blue") +
  geom_ribbon(aes(ymin = psiA_L, ymax = psiA_U), alpha = 0.2) +
  labs(title = NA,
       x = "Year", y = "Occupancy Index") +
  theme_minimal()


