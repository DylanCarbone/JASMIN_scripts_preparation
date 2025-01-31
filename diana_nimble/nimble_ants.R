library(tidyverse)
library(nimble)
library(MCMCvis)

task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))#39
#task.id = 4

set.seed(10)

# model directory --------------------------------------------------------------

setwd("./nimble_tests")

# load in the dataset ----------------------------------------------------------

data <- readRDS("monad_occupancy_dataset_ants.rds") %>%
                  rename(Species = tik,
                         SiteID = GRIDREF,
                         Date = lower_date)

# add date info ----------------------------------------------------------------

data <- data %>%
        mutate(Date = as.Date(Date),
              yday = lubridate::yday(Date),
              Year = lubridate::year(Date))

# site summary -----------------------------------------------------------------

siteSummary <- data %>%
                group_by(SiteID) %>%
                summarise(nuYears = length(unique(Year))) %>%
                ungroup() %>%
                filter(nuYears > 1)

data <- data %>%
          filter(SiteID %in% siteSummary$SiteID)

# species summary --------------------------------------------------------------

speciesSummary <- data %>%
                    group_by(Species) %>%
                    summarise(nuSiteID = length(unique(SiteID)), 
                              nuRecs = length(Species)) %>%
                    arrange(desc(nuRecs))

#choose species for run --------------------------------------------------------

#how many tasks do we have?
allSpecies <- sort(speciesSummary$Species[speciesSummary$nuRecs>50])
length(allSpecies)#46
myspecies <- allSpecies[task.id]
print(myspecies)

# define a visit ---------------------------------------------------------------

data <- data %>%
  mutate(visit = paste(Date, SiteID, sep="_"))

# organise species and visits --------------------------------------------------

#get occurence matrix  - detection and non-detection
occMatrix <- reshape2::acast(data, visit~Species, value.var="Species", fun=length)
occMatrix[occMatrix>0] <- 1

#get list length
visit_df <- data %>% 
  dplyr::group_by(visit, Year, SiteID) %>%
  dplyr::summarise(nuSpecies=length(unique(Species)),
                   nuRecords=length(Species)) %>%
  ungroup()

visit_df$yearIndex <- as.numeric(as.factor(visit_df$Year))
visit_df$siteIndex <- as.numeric(as.factor(visit_df$SiteID))

# order data
visit_df <- arrange(visit_df, visit)

# get data for a species -------------------------------------------------------

#check all aligns
all(row.names(occMatrix)==visit_df$visit)
visit_df$Species <- occMatrix[,myspecies]
table(visit_df$Species)

# model params -----------------------------------------------------------------

n.chains = 3
ni <- 32000 
nb <- 30000 
nt <- 5

# bundle data ------------------------------------------------------------------

nimbleConstants <- list(nsite = length(unique(visit_df$siteIndex)),
                        nyear = length(unique(visit_df$yearIndex)),
                        nvisit = nrow(visit_df),
                        site = visit_df$siteIndex,
                        year = visit_df$yearIndex,
                        LL = as.numeric(visit_df$nuSpecies))

#data
nimbleData <- list(y = as.numeric(visit_df$Species))

# initial values ---------------------------------------------------------------

#specify inits
z_inits <- reshape2::acast(visit_df, siteIndex ~ yearIndex, value.var="Species", fun = max)
#ignore warning message
z_inits[is.na(z_inits)] <- 0
z_inits[!z_inits %in% c(0,1)]<- 0

nimbleInits <- list(z = z_inits)

# nimble model -----------------------------------------------------------------

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


# fit model ------------------------------------------------------------------

#start <- print(Sys.time())
# mcmc.out <- nimbleMCMC(code = mynimbleCode,
#                        constants = nimbleConstants,
#                        data = nimbleData,
#                        inits = nimbleInits,
#                        nchains = n.chains,
#                        niter = ni, nburnin = nb, thin = nt,
#                        #samples = FALSE, summary = TRUE,
#                        samplesAsCodaMCMC = TRUE,
#                        monitors = c('sd.s','year.fixed','prop.p','beta.p'))
#end <- print(Sys.time())
#print(end-start)

# run line by line -------------------------------------------------------------

start <- print(Sys.time())

myModel <- nimbleModel(code = mynimbleCode,
                       data = nimbleData,
                       constants = nimbleConstants,
                       inits = nimbleInits)

CmyModel <- compileNimble(myModel)

myMCMC <- buildMCMC(CmyModel, monitors = c('sd.s','year.fixed','prop.p','beta.p'))

CmyMCMC <- compileNimble(myMCMC)

mid <- print(Sys.time())

results <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, 
                   nchains = n.chains)

end <- print(Sys.time())

#save output
temp <- MCMCsummary(object = results, round = 5)

#add times
temp$Stage1 <- mid-start
temp$Stage2 <- end-mid

saveRDS(temp, file=paste0("mcmc.summary_",myspecies,".rds"))


# # Processing outputs -----------------------------------------------------------
# outputs were downloaded and processing on local PC:
# out <- list.files("outputs", full.names = TRUE) %>%
#         map_dfr(readRDS)
# 
# #in seconds
# summary(as.numeric(out$Stage1)/(60*60))
# 
# #in hours(not sure why)
# summary(as.numeric(out$Stage2))

# end --------------------------------------------------------------------------
