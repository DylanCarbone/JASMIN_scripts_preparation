# setup the session and R ...
module load jaspy # For loading python workspace
module load jasr # For loading R workspace
R

library(parallel)
library(rslurm)
library(nimble)
library(tidyverse)
library(MCMCvis)
library(dplyr)

# Set working directory
setwd("dylcar_explore_occ_user")

# Read and preprocess data
data <- readRDS("ants/monad_occupancy_dataset_ants.rds") %>%
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

# Summarise species
speciesSummary <- data %>%
  group_by(Species) %>%
  summarise(nuSiteID = length(unique(SiteID)),
            nuRecs = length(Species)) %>%
  arrange(desc(nuRecs))

# Select species
allSpecies <- sort(speciesSummary$Species[speciesSummary$nuRecs > 50])

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
         siteIndex = as.numeric(as.factor(SiteID))) %>%
  arrange(visit)

# Model parameters
n.chains <- 3
ni <- 32000#32000
nb <- 30000#30000
nt <- 5

# Detect cores and calculate nodes
cores_for_clustering <- 15 #15
n_nodes_required <- ceiling(length(allSpecies) / cores_for_clustering)
species_to_node <- rep(1:n_nodes_required, each = cores_for_clustering, length.out = length(allSpecies))

nimbleConstants <- list(
  nsite = length(unique(visit_df$siteIndex)),
  nyear = length(unique(visit_df$yearIndex)),
  nvisit = nrow(visit_df),
  site = visit_df$siteIndex,
  year = visit_df$yearIndex,
  LL = as.numeric(visit_df$nuSpecies)
)

mynimbleCode <- nimbleCode({
  
  # State model
  for (i in 1:nsite){
    for (t in 1:nyear){
      
      # assume annual site occupancy can be reasonably explained by random effects for year and site
      logit(psi[i,t]) <-  year.fixed[t] + site.random[i]
      z[i,t] ~ dbern(psi[i,t])
      
    }
  }
  
  # Observation Model
  for(j in 1:nvisit) {
    
    y[j] ~ dbern(Py[j])
    y.sim[j] ~ dbern(Py[j])
    Py[j]<- z[site[j],year[j]]*p[j]
    
    #detection model
    #detection is affect by some covariates (list length)
    logit(p[j]) <- alpha.p + beta.p * LL[j]
    
  }
  
  # priors
  
  # State model
  
  #year effect
  for(t in 1:nyear){
    year.fixed[t] ~ dnorm(0, sd = 1.5) 
    logit(annual.occ[t]) <- year.fixed[t]
  }
  
  #site effect
  for(i in 1:nsite){
    site.random[i] ~ dnorm(0, sd = sd.s)
  }
  sd.s ~ T(dt(0, 1, 1), 0, 10)
  
  # Observation model
  
  #covs
  prop.p ~ dunif(0,1)
  logit(alpha.p) <- prop.p 
  
  beta.p ~ dnorm(0, 0.01)
  
})

# We'd need to create an initialisation dataset
first_species = allSpecies[1]
initial_df <- visit_df
initial_df$Species <- occMatrix[,first_species]

nimble_init_data <- list(y = as.numeric(initial_df$Species))

#specify inits
z_inits <- reshape2::acast(initial_df, siteIndex ~ yearIndex, value.var="Species", fun = max)
#ignore warning message
z_inits[is.na(z_inits)] <- 0
z_inits[!z_inits %in% c(0,1)]<- 0

# Main function
Initialise_NIMBLE_model_species <- function(node_id) {
  node_start_time <- Sys.time()

  # Initialisation
  nimbleInits <- list(z = z_inits)
  myModel <- nimbleModel(
    code = mynimbleCode,
    data = nimble_init_data,
    constants = nimbleConstants,
    inits = nimbleInits
  )

  print("1")
  
  
  CmyModel <- compileNimble(myModel)

  print("2")

  myMCMC <- buildMCMC(CmyModel, monitors = c('sd.s', 'year.fixed', 'prop.p', 'beta.p'))

  print("3")
  
  CmyMCMC <- compileNimble(myMCMC)

  print("4")

  intialisation_run_time <- as.numeric(difftime(Sys.time(), node_start_time, units = "mins"))

  print("5")

  # Species for the node
  node_species <- allSpecies[which(species_to_node == node_id)]

  print("6")
  
  logs <- parallel::mclapply(seq_along(node_species), function(core_i) {
    core_start_time <- Sys.time()

    species_name <- node_species[core_i]
    visit_df$Species <- occMatrix[, species_name]

    # Prepare data for this species
    nimbleData <- list(y = as.numeric(visit_df$Species))

    z_inits <- reshape2::acast(visit_df, siteIndex ~ yearIndex, value.var = "Species", fun = max)
    z_inits[is.na(z_inits)] <- 0
    z_inits[!z_inits %in% c(0, 1)] <- 0

    print("7")

    CmyModel$setData(nimbleData)
    CmyModel$setInits(list(z = z_inits))

    print("8")

    tryCatch({
      # Run the model
      results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = n.chains)
      temp <- MCMCsummary(object = results, round = 5)
      saveRDS(temp, file = paste0("mcmc.summary_", species_name, ".rds"))

      output_size <- as.numeric(object.size(temp)) / (1024^3) # Size in GB
    }, error = function(e) {
      output_size <- NA
      writeLines(paste0("Error for species: ", species_name, "\n", conditionMessage(e)), paste0(species_name, "_error.txt"))
    })

    core_end_time <- Sys.time()

    log_entry <- data.frame(
      taxa_group = "Ants",
      species_name = species_name,
      JASMIN = TRUE,
      queue = "par-single",
      n_nodes_requested = n_nodes_required,
      n_node_cores = cores_for_clustering,
      NIMBLE_initialisation_run_time = intialisation_run_time,
      node_iteration = node_id,
      node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
      node_end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      node_run_time = as.numeric(difftime(Sys.time(), node_start_time, units = "hours")),
      core_iteration = core_i,
      core_start_time = format(core_start_time, "%Y-%m-%d %H:%M:%S"),
      core_end_time = format(core_end_time, "%Y-%m-%d %H:%M:%S"),
      core_run_time = as.numeric(difftime(core_end_time, core_start_time, units = "hours")),
      output_object_size_GB = output_size,
      stringsAsFactors = FALSE
    )

    return(log_entry)
  }, mc.cores = cores_for_clustering)

  # Combine logs and save
  combined_logs <- do.call(rbind, logs)
  write.csv(combined_logs, paste0("node_", node_id, "_log.csv"))
}

# Generate the job name with the current date
jobname <- paste0('dylcar_explore_occ_run_NIMBLE_', format(Sys.Date(), "%d_%m_%Y"))

# Slurm job submission
sjob <- slurm_apply(
  f = Initialise_NIMBLE_model_species,
  params = data.frame(node_id = 1:n_nodes_required),
  jobname = jobname,
  nodes = n_nodes_required,
  cpus_per_node = cores_for_clustering,
  submit = TRUE,
  global_objects = c("z_inits", "allSpecies", "species_to_node", "nimbleConstants", "mynimbleCode", "nimble_init_data", "visit_df", "occMatrix", "ni", "nb", "nt", "n.chains", "cores_for_clustering"),
  slurm_options = list(time = "10:00:00", mem = 40000, partition = "par-single", error = "%a.err")
)

# Decide whether to make outputs match outputs of sparta OR change the BRCindicators inputs to match the NIMBLE outputs.
# Francesca's modified functions are probably based on the modifications that louise and Kath made.
# setup conversation with Diana and Nick to compare workflows.
# Send Kath objects and she may be able to take a lok
