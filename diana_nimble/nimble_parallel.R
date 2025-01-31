library(parallel)
library(nimble)
library(MCMCvis)

print(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# parLapply approach -----------------------------------------------------------

#option 1 from: https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html 

print(detectCores())

this_cluster <- makeCluster(4)

set.seed(10120)

useWAIC <- TRUE

# Simulate some data
myData <- rgamma(1000, shape = 0.4, rate = 0.8)

myCode <- nimbleCode({
  a ~ dunif(0, 100)
  b ~ dnorm(0, 100)
  
  for (i in 1:length_y) {
    y[i] ~ dgamma(shape = a, rate = b)
  }
})

# Create a function with all the needed code
run_MCMC_allcode <- function(seed, data, code, monitors) {
  
  library(nimble)
  
  myModel <- nimbleModel(code = code,
                         data = list(y = data),
                         constants = list(length_y = 1000),
                         inits = list(a = 0.5, b = 0.5))
  
  CmyModel <- compileNimble(myModel)
  
  myMCMC <- buildMCMC(CmyModel, monitors = monitors)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, niter = 10000, setSeed = seed)
  
  return(results)
  
}

start <- print(Sys.time())

chain_output <- parLapply(cl = this_cluster, 
                          X = 1:4, 
                          fun = run_MCMC_allcode, 
                          data = myData, code = myCode,
                          monitors = c("a", "b"))

#stopCluster(this_cluster)
end <- print(Sys.time())
print(end-start)

MCMCsummary(object = chain_output, round = 5)

# non-parallelized -------------------------------------------------------------

start <- print(Sys.time())

mcmc.out <- nimbleMCMC(code = myCode,
                       data = list(y = myData),
                       constants = list(length_y = 1000),
                       inits = list(a = 0.5, b = 0.5),
                       nchains = 4, 
                       niter = 10000, 
                       samplesAsCodaMCMC = TRUE,
                       monitors = c('a','b'))

end <- print(Sys.time())

print(end-start)

# # allowing different starting values per chain -------------------------------
# 
# this_cluster <- makeCluster(4) # Use only 2 for simplicity
# 
# set.seed(10120)
# 
# # Simulate some data
# 
# myData <- list(y = rgamma(1000, shape = 0.4, rate = 0.8))
# myConstants <- list(length_y = 1000)
# 
# myCode <- nimbleCode({
#   a ~ dunif(0, 100)
#   b ~ dnorm(0, 100)
#   
#   for (i in 1:length_y) {
#     y[i] ~ dgamma(shape = a, rate = b)
#   }
# })
# 
# # Option 1: make global lists and use the X input from parLapply
# # to look up elements.
# 
# run_MCMC_allcode <- function(X, data, code, constants) {
#   
#   library(nimble)
#   
#   inputs <- inputsList[[X]]
#   inits <- inputs$inits
#   seed <- inputs$seed
#   
#   myModel <- nimbleModel(code = code,
#                          data = data,
#                          constants = constants,
#                          inits = inits)
#   
#   
#   CmyModel <- compileNimble(myModel)
#   
#   myMCMC <- buildMCMC(CmyModel)
#   
#   CmyMCMC <- compileNimble(myMCMC)
#   
#   results <- runMCMC(CmyMCMC, niter = 10000, setSeed = seed)
#   
#   return(results)
#   
# }
# 
# 
# inputsList <- list(
#   list(seed = 123, # The seeds are arbitrary numbers
#        inits = list(a = 0.5, b = 0.5)),
#   list(seed = 234,
#        inits = list(a = 0.3, b = 0.3)),
#   list(seed = 345,
#        inits = list(a = 0.3, b = 0.8)),
#   list(seed = 456,
#        inits = list(a = 0.3, b = 0.2)))
# 
# clusterExport(this_cluster, c("inputsList"))
# 
# chain_output <- parLapply(cl = this_cluster, X = 1:4,
#                           fun = run_MCMC_allcode,
#                           data = myData,
#                           code = myCode,
#                           constants = myConstants)
# 
# MCMCsummary(object = chain_output, round = 5)
# 
# # pbmcapply approach -----------------------------------------------------------
# 
# myData <- rgamma(1000, shape = 0.4, rate = 0.8)
# 
# myCode <- nimbleCode({
#   a ~ dunif(0, 100)
#   b ~ dnorm(0, 100)
#   
#   for (i in 1:length_y) {
#     y[i] ~ dgamma(shape = a, rate = b)
#   }
# })
# 
# model <- nimbleModel(code = myCode,
#                      constants = list(length_y = 1000),
#                      data = list(y = myData),
#                      inits = list(a = 0.1, b = 0.2))
# 
# occMCMC <- buildMCMC(model,
#                      monitors = c("a","b"),
#                      thin = 1,
#                      useConjugacy = FALSE)
# 
# Cmodel <- compileNimble(model)
# 
# CoccMCMC <- compileNimble(occMCMC, project = model)
# 
# av_cores <- parallel::detectCores() - 1
# 
# n.chain = 4
# 
# runMCMC_samples <- pbmcapply::pbmclapply(1:n.chain, function(i)
#   
#   runMCMC(
#     mcmc = CoccMCMC,
#     nburnin = n.burn,
#     niter = n.iter,
#     nchains = 1, 
#     samplesAsCodaMCMC = T),
#     mc.cores = av_cores)
# 
# runMCMC_samples