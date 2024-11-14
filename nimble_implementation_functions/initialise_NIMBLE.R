initialise_NIMBLE <- function() {

  dataConstants <- formattedData$dataConstants
  obsData <- formattedData$obsData

  # Example code to measure
  start_time <- Sys.time()  # Start time

  # kludge required: if maxSp is set to 1 (or zero) then problems arise later.
  if(maxSp < 2) maxSp <- 2

  # The number of species to be modelled
  nSpMod <- min(dim(obsData$y)[1], maxSp)

  if(is.null(n.burn)) n.burn <- round(n.iter / 2)

  # truncate the dataset if there are too many species
  if(dim(obsData$y)[1] > nSpMod){
    obsData$y <- obsData$y[1:nSpMod,]
    dataSumm$occMatrix <- dataSumm$occMatrix[1:nSpMod,,]
    dataSumm$stats <- dataSumm$stats[1:nSpMod,]
    dataConstants$nsp <- nSpMod
    print(paste('Warning: only the first', nSpMod, 'species will be used in modelling: others will be ignored'))
  }

  # put the species names in a vector
  spNames <- dimnames(formattedData$obsData$y)[[1]][1:nSpMod]

  # define the model parameters and which should be monitored
  modPars <- c('lam.0', 'psi.fs', "sd.psi", "alpha.0", "mu.alpha", "sd.alpha")
  
  if(inclPhenology) modPars <- c(modPars, "alpha.1", "beta1", "beta2")
  if(inclStateRE) modPars <- c(modPars, "sd.eta")
  
  if(!is.null(ListLen)){
    modPars <- c(modPars, "gamma.1")
    if(ListLen == "cat") modPars <- c(modPars, "gamma.2") 
  }

  if(all(is.logical(all.Pars))){
    if(all.Pars == TRUE) {
      params <- modPars
    } else {
      params <- c("lam.0")
    }
  } else {
    params <- all.Pars
    if(!all(params %in% modPars)){
      badPars <- setdiff(params, modPars)
      warning(paste(badPars, "not recognised"))
      if(all(params %in% badPars)) {
        stop("no valid parameters listed")
      } else {
        params <- setdiff(params, badPars)
      }
    }
  }
  print(paste("Monitoring:", params))

  params2 <- intersect(params, c("psi.fs", "alpha.0"))
  if(any(c("psi.fs", "alpha.0") %in% params)){
    params <- setdiff(params, c("psi.fs", "alpha.0"))
  }

  # step 1 define the model code
  modelcode <- defineModel_SS(inclPhenology = inclPhenology,
                              ListLen = ListLen,
                              inclStateRE = inclStateRE)

  init.vals <- list(z = dataSumm$occMatrix[1,,], 
                    lam.0 = logit(dataSumm$stats$naiveOcc[1] * 0.99), 
                    sd.psi = 0.5,
                    psi = rnorm(n=dataConstants$nyear, sd = 0.02),
                    mu.alpha = -1,
                    sd.alpha = 2
                   )

  if(inclPhenology){
    init.vals$beta1 <- 180
    init.vals$beta2 <- 50
    init.vals$alpha.1 <- 1
    init.vals$alpha.0 <- rnorm(n=dataConstants$nyear, mean=0, sd=2)
  } else {
    init.vals$alpha.0 <- rnorm(n=dataConstants$nyear, mean= -2, sd=2)
  }
  
  if(inclStateRE){
    init.vals$sd.eta <- 2
    init.vals$eta <- rnorm(n=dataConstants$nsite, mean=0, sd=2)
  }
  
  if(!is.null(ListLen)) {
    if(ListLen == "cont"){
      init.vals$gamma.1 <- 0.1
    } else if(ListLen == "cat"){
      init.vals$gamma.1 <- 0.1
      init.vals$gamma.2 <- 0.1
    } else {
      stop("invalid List Length option")
    }
  }

  # step 2 create an operational model from NIMBLE/BUGS code
  model <- nimbleModel(code = modelcode,
                       constants = dataConstants[!names(dataConstants) %in% "nsp"],
                       data = lapply(obsData, function(x) x[1,]),
                       inits = init.vals)

  print("Nimble initialisation complete")

  # step 3 build an MCMC object using buildMCMC()
  occMCMC <- buildMCMC(model,
                       monitors = params,
                       monitors2 = params2,
                       useConjugacy = FALSE)
  print("Model MCMC build complete")

  # step 4 compile the model
  Cmodel <- compileNimble(model)
  print("Model compilation complete")

  # compile the MCMC object
  CoccMCMC <- compileNimble(occMCMC, project = model)
  print("MCMC compilation complete")

  # Record execution time
  end_time <- Sys.time()    
  execution_time <- end_time - start_time

  transferable_outputs <- list(nSpMod = nSpMod, data = obsData, dataSumm = dataSumm, n.iter = n.iter,
   n.burn = n.burn, Cmodel = Cmodel, CoccMCMC = CoccMCMC, mon2 = length(params2)>0)

  saveRDS(transferable_outputs, file = "transferable_outputs.rds")

  print(paste("Execution time:", execution_time))
  write(paste("Execution time:", execution_time), file = "execution_time.txt")
  
}