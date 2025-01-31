library(tidyverse)
library(nimble)
library(igraph)
library(coda)
library(MCMCvis)
library(parallel)

#set model parameters
#n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) 
n.cores = 3
#ni <- 20  ;   nb <- 10   ;   nt <- 2 
#ni <- 10000  ;   nb <- 5000   ;   nt <- 20 
ni <- 200000 
nb <- ni*(3/4)   
nt <- 50

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
         siteIndex = as.numeric(as.factor(SiteID)))

nimbleConstants <- list(nsite = length(unique(visit_df$siteIndex)),
                        nyear = length(unique(visit_df$yearIndex)),
                        nvisit = nrow(visit_df),
                        site = visit_df$siteIndex,
                        year = visit_df$yearIndex,
                        m_p = m_p, # number of parameters
                        Covs_p = Covs_p, #design matrix for detection
                        X_spline = X_spline,
                        X_pred = X_pred,
                        predIndices = as.numeric(as.factor(predIndices)),
                        npredIndices = length(predIndices),
                        S1 = S1,
                        zero = zero)

myspecies <- allSpecies[species_i]
print(myspecies)

visit_df$Species <- occMatrix[,myspecies]
nimbleData <- list(y = as.numeric(visit_df$Species))

#data
nimbleData <- list(y = as.numeric(visit_df$Species))

#specify inits
z_inits <- reshape2::acast(visit_df, siteIndex ~ yearIndex, value.var="Species", fun = max)
#ignore warning message
z_inits[is.na(z_inits)] <- 0
z_inits[!z_inits %in% c(0,1)]<- 0

nimbleInits <- function() {
  list(z = z_inits,
       rho = rnorm(length(init_rho), 0, 2),
       b = init_b,
       sd.s = runif(1, 0, 2),
       beta_p = rnorm(m_p, 0, 1),
       site.random = rnorm(n = nimbleConstants$nsite, 0, 1),
       y.sim = nimbleData$y)
  }


mynimbleCode <- nimbleCode({

  # State/Ecological model
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
    #detection is affect by some covariates (list length) and day of year(modelled with a spline)
    logit(p[j]) <-  inprod(beta_p[1:m_p], Covs_p[j,1:m_p]) + inprod(b[1:8], X_spline[j,1:8])

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
  beta_p[1] ~ dnorm(0, sd = 1.5) #overall intercept
  
  for(i in 2:m_p){
    beta_p[i] ~ dnorm(0, sd = 10)
  }

  ## prior for s(yday)
  K1[1:8,1:8] <- S1[1:8,1:8] * lambda[1] 
  b[1:8] ~ dmnorm(zero[1:8],K1[1:8,1:8])

  ## smoothing parameter priors
  for (i in 1:1) {
    rho[i] ~ T(dnorm(0, sd = 2),-8, 8)
    lambda[i] <- exp(rho[i])
  }

  #predict the flight curve
  for(i in 1:npredIndices){
    logit(occ.day[predIndices[i]]) <- beta_p[1] + inprod(b[1:8], X_pred[i,1:8])
  }
  
  #predict annual occupancy change
  for(i in 1:(nyear-1)){
    occ.change[i] <- annual.occ[i+1]/annual.occ[i] 
  }
  
})

run_MCMC_allcode <- function(X, data, constants, code, monitors, ni, nb, nt) {
  
  inputs <- inputsList[[X]]
  inits <- inputs$inits
  seed <- inputs$seed
  
  myModel <- nimbleModel(code = code,
                         data = data,
                         constants = constants,
                         inits = inits)
  
  CmyModel <- compileNimble(myModel)
  
  myMCMC <- buildMCMC(CmyModel, monitors = monitors)
  
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, setSeed = seed)
  
  return(results)
  
}

inputsList <- list(
  list(seed = 123, inits = nimbleInits()),
  list(seed = 234, inits = nimbleInits()),
  list(seed = 345, inits = nimbleInits()))

this_cluster <- makeCluster(n.cores)
  
clusterExport(this_cluster, c("inputsList"))

print(detectCores())

start <- print(Sys.time())

chain_output <- parLapply(cl = this_cluster, 
                          X = 1:n.cores,
                          fun = run_MCMC_allcode,
                          data = nimbleData,
                          code = mynimbleCode,
                          constants = nimbleConstants,
                          monitors = c('rho','beta_p','sd.s',
                                       'occ.day','annual.occ','occ.change',
                                       'e.count','sim.count'),
                          ni = ni,
                          nb = nb,
                          nt = nt)


end <- print(Sys.time())

print(end-start)

#save output
temp <- MCMCsummary(object = chain_output, round = 5)
temp$Param <- row.names(temp)

saveRDS(temp, file=paste0("all_mcmc.summary_",myspecies,".rds"))