# This demonstrates how to “save” the model state, and internal state of an MCMC algorithm to disk.
# 
# Then, R can be restarted, and this saved “state” can be reloaded into a new (structurally identical) model and MCMC algorithm, allowing the original sampling chain to be “resumed” where it left off.
# 
# Function Definitions


getStateVariableNames <- function(samplerDef) {
  resetMethod <- body(samplerDef$reset)
  stateVars <- character()
  if(resetMethod[[1]] != '{') stop('something wrong')
  numLines <- length(resetMethod)
  for(i in 1:numLines) {
    if(i == 1) next
    thisLine <- resetMethod[[i]]
    if(thisLine[[1]] == '<<-') {
      LHS <- thisLine[[2]]
      if(!is.name(LHS)) stop('haven\'t dealt with non-name-LHS case yet')
      stateVars <- c(stateVars, as.character(LHS))
    }
    if('my_calcAdaptationFactor' %in% all.names(thisLine)) {
      stateVars <- c(stateVars, 'my_calcAdaptationFactor')
    }
  }
  setupMethod <- body(samplerDef$setup)
  if('empirSamp' %in% all.names(setupMethod)) stateVars <- c(stateVars, 'empirSamp')
  return(stateVars)
}

getModelState <- function(model) {
  modelVarNames <- model$getVarNames()
  modelVarValuesList <- vector('list', length(modelVarNames))
  names(modelVarValuesList) <- modelVarNames
  for(var in modelVarNames) {
    modelVarValuesList[[var]] <- model[[var]]
  }
  return(modelVarValuesList)
}

getMCMCstate <- function(conf, mcmc) {
  stateVarNamesList <- vector('list', length(conf$samplerConfs))
  mcmcStateValuesList <- vector('list', length(conf$samplerConfs))
  for(i in seq_along(conf$samplerConfs)) {
    samplerDef <- conf$getSamplerDefinition(i)
    theseStateNames <- getStateVariableNames(samplerDef)
    theseStateValuesList <- vector('list', length(theseStateNames))
    names(theseStateValuesList) <- theseStateNames
    for(j in seq_along(theseStateNames)) {
      if(is.nf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$timesAdapted,
                                            gamma1 = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$gamma1)
        } else
          theseStateValuesList[[j]] <- mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]
      }
      if(is.Cnf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted'),
                                            gamma1 = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1'))
        } else
          theseStateValuesList[[j]] <- valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], theseStateNames[j])
      }
    }
    mcmcStateValuesList[[i]] <- theseStateValuesList
  }
  return(mcmcStateValuesList)
}

getWAICstate <- function(mcmc) {
  waicStateNames <- c('delta1pWAICmat', 'delta2pWAICmat', 'finalized', 'logProbMat', 'lppdCurSumMat', 'lppdSumMaxMat',
                      'mcmcIter', 'meanpWAICmat', 'sspWAICmat')
  waicStateList <- vector('list', length(waicStateNames))
  names(waicStateList) <- waicStateNames
  for(nm in waicStateNames) {
    if(is.nf(mcmc)) {
      waicStateList[[nm]] <- mcmc$waicFun[[1]][[nm]]
    }
    if(is.Cnf(mcmc)) {
      waicStateList[[nm]] <- valueInCompiledNimbleFunction(mcmc$waicFun[[1]], nm)
    }
  }
  return(waicStateList)
}

setModelState <- function(model, modelState) {
  modelVarNames <- model$getVarNames()
  if(!identical(sort(modelVarNames), sort(names(modelState)))) stop('saved model variables don\'t agree')
  for(var in modelVarNames) {
    model[[var]] <- modelState[[var]]
  }
  invisible(model$calculate())
}

setMCMCstate <- function(conf, mcmc, mcmcState) {
  if(length(mcmcState) != length(conf$samplerConfs)) stop('saved mcmc samplers don\'t agree')
  for(i in seq_along(conf$samplerConfs)) {
    theseStateValuesList <- mcmcState[[i]]
    for(j in seq_along(theseStateValuesList)) {
      samplerStateName <- names(theseStateValuesList)[j]
      if(is.nf(mcmc)) {
        if(samplerStateName == 'my_calcAdaptationFactor') {
          mcmc$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$timesAdapted <- theseStateValuesList[[j]]$timesAdapted
          mcmc$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$gamma1 <- theseStateValuesList[[j]]$gamma1
        } else {
          mcmc$samplerFunctions$contentsList[[i]][[samplerStateName]] <- theseStateValuesList[[samplerStateName]]
        }
      }
      if(is.Cnf(mcmc)) {
        if(samplerStateName == 'my_calcAdaptationFactor') {
          invisible(valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted', theseStateValuesList[[j]]$timesAdapted))
          invisible(valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1', theseStateValuesList[[j]]$gamma1))
        } else {
          invisible(valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], samplerStateName, theseStateValuesList[[samplerStateName]]))
        }
      }
    }
  }
}

setWAICstate <- function(mcmc, waicState) {
  for(nm in names(waicState)) {
    if(is.nf(mcmc)) {
      mcmc$waicFun[[1]][[nm]] <- waicState[[nm]]
    }                       
    if(is.Cnf(mcmc)) {
      invisible(valueInCompiledNimbleFunction(mcmc$waicFun[[1]], nm, waicState[[nm]]))
    }
  }
}


# Usage Example
library(nimble)

# Define a NIMBLE model object, and an MCMC algorithm to operate on it. Optionally, choose to calculate WAIC.
# 
# The saving / loading state will work for most any NIMBLE MCMC algorithms, but if you have questions about your use case, then please email Daniel at dbt1@williams.edu.



useWAIC <- TRUE

code <- nimbleCode({
  for(i in 1:N) {
    x[i] ~ dnorm(0, 1)
  }
  b ~ dbern(0.5)
  c ~ dnorm(0, 1)
  y ~ dnorm(c, 1)
  p ~ dnorm(0, 1)
})

N <- 8
constants <- list(N = N)
data <- list(y = 0)
inits <- list(x = rep(0, N),
              b = 0,
              c = 0,
              p = 0)

## create model object
Rmodel <- nimbleModel(code, constants, data, inits)
## Defining model
## Building model
## Setting data and initial values
## Running calculate on model
##   [Note] Any error reports that follow may simply reflect missing values in model variables.
## Checking model sizes and dimensions
Rmodel$calculate()
## [1] -10.80147
## configure & build MCMC
conf <- configureMCMC(Rmodel, enableWAIC = useWAIC)
## ===== Monitors =====
## thin = 1: b, c, p, x
## ===== Samplers =====
## posterior_predictive sampler (10)
##   - x[]  (8 elements)
##   - b
##   - p
## conjugate sampler (1)
##   - c
conf$removeSamplers(c('b', 'x'))
conf$addSampler('b', 'binary')
conf$addSampler('x[1]', 'RW')
conf$addSampler('x[2]', 'slice')
conf$addSampler('x[3:5]', 'RW_block')
##   [Note] Assigning an RW_block sampler to nodes with very different scales can result in low MCMC efficiency.  If all nodes assigned to RW_block are not on a similar scale, we recommend providing an informed value for the "propCov" control list argument, or using the AFSS sampler instead.
conf$addSampler('x[6:8]', 'AF_slice')
conf$printSamplers()
## [1] conjugate_dnorm_dnorm_identity sampler: c
## [2] posterior_predictive sampler: p
## [3] binary sampler: b
## [4] RW sampler: x[1]
## [5] slice sampler: x[2]
## [6] RW_block sampler: x[3:5]
## [7] AF_slice sampler: x[6:8]
Rmcmc <- buildMCMC(conf)
##   [Note] Reordering posterior_predictive samplers to execute last
## compile
Cmodel <- compileNimble(Rmodel)
## Compiling
##   [Note] This may take a minute.
##   [Note] Use 'showCompilerOutput = TRUE' to see C++ compilation details.
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
## Compiling
##   [Note] This may take a minute.
##   [Note] Use 'showCompilerOutput = TRUE' to see C++ compilation details.
## execute MCMC for 10000 iterations
set.seed(0)
Cmcmc$run(10000)
## |-------------|-------------|-------------|-------------|
## |-------------------------------------------------------|
## NULL


# Now, we’ll “extract” and save the model state, and the internal state of the MCMC.
# 
# For saving the MCMC state, you’ll need both the compiled MCMC object (Cmcmc), and also the MCMC configuration object (conf).

stateList <- list(modelState = getModelState(Cmodel),
                  mcmcState = getMCMCstate(conf, Cmcmc))
if(useWAIC)
  stateList[['waicState']] <- getWAICstate(Cmcmc)

fileName <- "first_run.RDS"
saveRDS(stateList, file = fileName)

# Here, you could restart R.
# 
# Now we’ll make new (compiled) model and MCMC objects, then load this previous “state” back into them.


library(nimble)
## construct an *identical* model object
Rmodel2 <- nimbleModel(code, constants, data, inits)
## Defining model
## Building model
## Setting data and initial values
## Running calculate on model
##   [Note] Any error reports that follow may simply reflect missing values in model variables.
## Checking model sizes and dimensions
## construct an *identical* MCMC algorithm
conf2 <- configureMCMC(Rmodel2, enableWAIC = useWAIC)
## ===== Monitors =====
## thin = 1: b, c, p, x
## ===== Samplers =====
## posterior_predictive sampler (10)
##   - x[]  (8 elements)
##   - b
##   - p
## conjugate sampler (1)
##   - c
conf2$removeSamplers(c('b', 'x'))
conf2$addSampler('b', 'binary')
conf2$addSampler('x[1]', 'RW')
conf2$addSampler('x[2]', 'slice')
conf2$addSampler('x[3:5]', 'RW_block')
##   [Note] Assigning an RW_block sampler to nodes with very different scales can result in low MCMC efficiency.  If all nodes assigned to RW_block are not on a similar scale, we recommend providing an informed value for the "propCov" control list argument, or using the AFSS sampler instead.
conf2$addSampler('x[6:8]', 'AF_slice')
conf2$printSamplers()
## [1] conjugate_dnorm_dnorm_identity sampler: c
## [2] posterior_predictive sampler: p
## [3] binary sampler: b
## [4] RW sampler: x[1]
## [5] slice sampler: x[2]
## [6] RW_block sampler: x[3:5]
## [7] AF_slice sampler: x[6:8]
Rmcmc2 <- buildMCMC(conf2)
##   [Note] Reordering posterior_predictive samplers to execute last
## compile
Cmodel2 <- compileNimble(Rmodel2)
## Compiling
##   [Note] This may take a minute.
##   [Note] Use 'showCompilerOutput = TRUE' to see C++ compilation details.
Cmcmc2 <- compileNimble(Rmcmc2, project = Rmodel2)
## Compiling
##   [Note] This may take a minute.
##   [Note] Use 'showCompilerOutput = TRUE' to see C++ compilation details.


# Now, we’ll “restore” the previous state into the new model (Cmodel2) and the new MCMC (Cmcmc2).


## load the saved "state" file
stateList <- readRDS(file = fileName)

modelState <- stateList$modelState
mcmcState <- stateList$mcmcState
if(useWAIC)
  waicState <- stateList$waicState

## restore the saved "state" into the new model and new MCMC
setModelState(Cmodel2, modelState)
setMCMCstate(conf2, Cmcmc2, mcmcState)
if(useWAIC)
  setWAICstate(Cmcmc2, waicState)


# Finally, now you can “resume” the MCMC run, exactly where it left off.
# 
# Important: Make certain to use the argument reset = FALSE, so the MCMC state is not reset to its initial settings for a new, fresh, MCMC run. If calculating WAIC, also set resetWAIC = FALSE.
# 
# This option, reset = FALSE, is only available through the mcmc$run(...) invocation of the MCMC algorithm.



## continue MCMC run
Cmcmc2$run(10000, reset = FALSE, resetWAIC = FALSE)
## |-------------|-------------|-------------|-------------|
## |-------------------------------------------------------|
## NULL


That’s it! You did it!
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  