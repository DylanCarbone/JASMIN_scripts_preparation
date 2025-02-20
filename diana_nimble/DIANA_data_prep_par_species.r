# Load necessary modules
module load jaspy  # For loading Python workspace ??
module load jasr    # For loading R workspace ??

R

library(parallel)
library(nimble)
library(tidyverse)
library(MCMCvis)
library(dplyr)
library(Matrix)

# Model parameters
n.chains <- 3
ni <- 100 # 32000
nb <- 90 # 30000
nt <- 5
min.Recs <- 50 # 10 # number of records per species for inclusion

# Number of cores and number of species per core
cores_for_clustering <- 5
species_to_core <- 3

# Set working directory
setwd("dylcar_explore_occ_user")

# Read and preprocess data
data <- readRDS("monad_occupancy_dataset_ants.rds") %>% 
  rename(Species = tik, SiteID = GRIDREF, Date = lower_date) %>%
  mutate(Date = as.Date(Date),
         yday = lubridate::yday(Date),
         Year = lubridate::year(Date))

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

n_nodes_required <- ceiling(length(allSpecies) / (cores_for_clustering * species_to_core))
species_to_node <- rep(1:n_nodes_required, each = (cores_for_clustering * species_to_core), length.out = length(allSpecies))

# Prepare data
data <- data %>%
  mutate(visit = paste(Date, SiteID, sep = "_"))

occMatrix <- reshape2::acast(data, visit ~ Species, value.var = "Species", fun = length)
occMatrix[occMatrix > 0] <- 1

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
  for (i in 1:nsite){
    for (t in 1:nyear){
      logit(psi[i,t]) <-  year.fixed[t] + site.random[i]
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
jobname <- paste0('dylcar_explore_occ_run_DIANA_NIMBLE_ANTS_par_species', "_", format(Sys.Date(), "%d_%m_%Y"))
dir.create(jobname)

saveRDS(allSpecies, file.path(jobname, "allSpecies.RDS"))
saveRDS(nimbleConstants, file.path(jobname, "nimbleConstants.RDS"))
saveRDS(mynimbleCode, file.path(jobname, "mynimbleCode.RDS"))
saveRDS(visit_df, file.path(jobname, "visit_df.RDS"))
saveRDS(occMatrix, file.path(jobname, "occMatrix.RDS"))
saveRDS(ni, file.path(jobname, "ni.RDS"))
saveRDS(nb, file.path(jobname, "nb.RDS"))
saveRDS(nt, file.path(jobname, "nt.RDS"))
saveRDS(n.chains, file.path(jobname, "n.chains.RDS"))
saveRDS(species_to_node, file.path(jobname, "species_to_node.RDS"))

# setwd(jobname)

# Please request:

#array size
paste0("0-", length(allSpecies)-1)

# Cores per node
cores_for_clustering + 1
