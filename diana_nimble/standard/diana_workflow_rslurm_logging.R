# setup the session and R ...
module load jaspy # For loading python workspace ??
module load jasr # For loading R workspace ??

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
min.Recs <- 50 # 50 # number of records per species for inclusion

# Set working directory
setwd("dylcar_explore_occ_user")

# load("butterfly_data/BNM_2022.rdata")
# butterfly_data = taxa_data %>%
# rename(tik = CONCEPT, lower_date = TO_STARTDATE, GRIDREF = TO_GRIDREF)

# saveRDS(butterfly_data, "formatted_butterfly_data.rds")  

# Read and preprocess data
data <- readRDS("formatted_butterfly_data.rds") %>% ### HERE
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
allSpecies <- sort(speciesSummary$Species[speciesSummary$nuRecs > min.Recs])

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

NIMBLE_run = function(species_i){

  node_start_time = Sys.time()

  myspecies <- allSpecies[species_i]
  print(myspecies)

  visit_df$Species <- occMatrix[,myspecies]
  nimbleData <- list(y = as.numeric(visit_df$Species))

  #specify inits
  z_inits <- reshape2::acast(visit_df, siteIndex ~ yearIndex, value.var="Species", fun = max)

  #ignore warning message
  z_inits[is.na(z_inits)] <- 0
  z_inits[!z_inits %in% c(0,1)]<- 0

  nimbleInits <- list(z = z_inits)

  myModel <- nimbleModel(code = mynimbleCode,
                        data = nimbleData,
                        constants = nimbleConstants,
                        inits = nimbleInits)

  CmyModel <- compileNimble(myModel)

  myMCMC <- buildMCMC(CmyModel, monitors = c('sd.s','year.fixed','prop.p','beta.p'))

  CmyMCMC <- compileNimble(myMCMC)

  intialisation_run_time <- as.numeric(difftime(Sys.time(), node_start_time, units = "hours"))

  results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, 
                    nchains = n.chains)

  #save output
  temp <- MCMCsummary(object = results, round = 5)

  saveRDS(temp, file=paste0("mcmc.summary_",myspecies,".rds"))

  log_entry <- data.frame(
    taxa_group = "Butterflies", ### HERE
    species_name = myspecies,
    JASMIN = TRUE,
    queue = "short-serial",
    n_nodes_requested = length(allSpecies),
    NIMBLE_initialisation_run_time = intialisation_run_time,
    node_start_time = format(node_start_time, "%Y-%m-%d %H:%M:%S"),
    node_end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    node_run_time = as.numeric(difftime(Sys.time(), node_start_time, units = "hours")),
    stringsAsFactors = FALSE
  )

  write.csv(log_entry, paste0(myspecies, "_log.csv"))

}

# Generate the job name with the current date
jobname <- paste0('dylcar_explore_occ_run_DIANA_NIMBLE_BUTTERFLIES', "_", format(Sys.Date(), "%d_%m_%Y")) ### HERE

# NOTE THE CURRENT DATASET LOADED IS ANTS
# Slurm job submission
sjob <- slurm_apply(
  f = NIMBLE_run,
  params = data.frame(species_i = 1:length(allSpecies)),
  jobname = jobname,
  nodes = length(allSpecies),
  cpus_per_node = 1,
  submit = TRUE,
  global_objects = c("allSpecies", "nimbleConstants", "mynimbleCode", "visit_df", "occMatrix", "ni", "nb", "nt", "n.chains"),
  slurm_options = list(time = "168:00:00", mem = 60000, partition = "long-serial", error = "%a.err") ### HERE
)
