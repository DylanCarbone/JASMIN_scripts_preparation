# Load necessary modules
module load jaspy  # For loading Python workspace ??
module load jasr    # For loading R workspace ??

R

library(parallel)
library(rslurm)
library(nimble)
library(tidyverse)
library(MCMCvis)
library(dplyr)
library(Matrix)

# Model parameters
n.chains <- 3
ni <- 32000 # 32000
nb <- 30000 # 30000
nt <- 5
min.Recs <- 50 # 50 number of records per species for inclusion

# FOR TESTING ONLY
# Take a sample of rows. If NULL, it is ignored
sample_records = NULL

# Set working directory
setwd("dylcar_explore_occ_user")

# Read and preprocess data
data <- readRDS("formatted_butterfly_data.rds") %>% 
  rename(Species = tik, SiteID = GRIDREF, Date = lower_date) %>%
  mutate(Date = as.Date(Date),
         yday = lubridate::yday(Date),
         Year = lubridate::year(Date))

if (!is.null(sample_records)){
  data = data %>% sample_n(sample_records)
}

# Filter data for sites with records across multiple years
siteSummary <- data %>%
  group_by(SiteID) %>%
  summarise(nuYears = length(unique(Year))) %>%
  ungroup() %>%
  filter(nuYears > 1)

data <- data %>% filter(SiteID %in% siteSummary$SiteID)

# Summarize species
speciesSummary <- data %>%
  group_by(Species) %>%
  summarise(nuSiteID = length(unique(SiteID)),
            nuRecs = length(Species)) %>%
  arrange(desc(nuRecs))

# Select species
allSpecies <- sort(speciesSummary$Species[speciesSummary$nuRecs > min.Recs]) # 50

# Prepare data
data <- data %>%
  mutate(visit = paste(Date, SiteID, sep = "_"))

species_counts <- split(data, data$Species) %>%
  lapply(function(df){ table(df$visit)})

visit_df <- data %>%
  group_by(visit, Year, SiteID) %>%
  summarise(nuSpecies = length(unique(Species)),
            nuRecords = length(Species)) %>%
  ungroup()

visit_df <- visit_df %>%
  mutate(yearIndex = as.numeric(as.factor(Year)),
         siteIndex = as.numeric(as.factor(SiteID)))

nimbleConstants <- list(
  nsite = length(unique(visit_df$siteIndex)),
  nyear = length(unique(visit_df$yearIndex)),
  nvisit = nrow(visit_df),
  site = visit_df$siteIndex,
  year = visit_df$yearIndex,
  LL = as.numeric(visit_df$nuSpecies)
)

mynimbleCode <- nimbleCode({
  for (t in 1:nyear) {  # Outer loop over years

    psi.fs[t] <- sum(z[1:nsite, t]) / nsite  # Now defined once per year

    for (i in 1:nsite) {  # Inner loop over sites
      logit(psi[i,t]) <- year.fixed[t] + site.random[i]
      z[i,t] ~ dbern(psi[i,t])
    }
  }

  for(j in 1:nvisit) {
    y[j] ~ dbern(Py[j])
    y.sim[j] ~ dbern(Py[j])
    Py[j] <- z[site[j],year[j]] * p[j]
    logit(p[j]) <- alpha.p + beta.p * LL[j]
  }

  for(t in 1:nyear){
    year.fixed[t] ~ dnorm(0, sd = 1.5) 
    logit(annual.occ[t]) <- year.fixed[t]
  }

  for(i in 1:nsite){
    site.random[i] ~ dnorm(0, sd = sd.s)
  }
  sd.s ~ T(dt(0, 1, 1), 0, 10)

  prop.p ~ dunif(0,1)
  logit(alpha.p) <- prop.p 
  beta.p ~ dnorm(0, 0.01)
})

# Generate the job name with the current date
jobname <- paste0('dylcar_explore_occ_run_DIANA_NIMBLE_BUTTERFLIES_par_threads', "_", format(Sys.Date(), "%d_%m_%Y"))
jobname_full <- paste0("_rslurm_", jobname)

dir.create(jobname_full)

# Save each of the species counts as a seperate dataframe
for (species in names(species_counts)){

  dir.create(file.path(jobname_full, "species_data"), showWarnings = FALSE)

  saveRDS(species_counts[[species]], file.path(jobname_full, "species_data", paste0(species, ".RDS")))
}

saveRDS(allSpecies, file.path(jobname_full, "allSpecies.RDS"))
saveRDS(nimbleConstants, file.path(jobname_full, "nimbleConstants.RDS"))
saveRDS(mynimbleCode, file.path(jobname_full, "mynimbleCode.RDS"))
saveRDS(visit_df, file.path(jobname_full, "visit_df.RDS"))
saveRDS(ni, file.path(jobname_full, "ni.RDS"))
saveRDS(nb, file.path(jobname_full, "nb.RDS"))
saveRDS(nt, file.path(jobname_full, "nt.RDS"))
saveRDS(n.chains, file.path(jobname_full, "n.chains.RDS"))

# Please request:

#array size
paste0("0-", length(allSpecies)-1)

# cpus-per-task
n.chains + 1