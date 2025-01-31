single_species_model <- function(sp, spDat, dataSumm,
                                       n.iter, n.burn, n.thin, n.chain,
                                       Cmodel, CoccMCMC, mon2 = FALSE){

        # apparent occupancy for this species
        Z <- dataSumm$occMatrix[sp,,]

        # write an informative message about this species' data
        nS <- sum(rowSums(Z, na.rm=T)>0)
        spName <- dataSumm$stats$species[sp]
        print(paste0("Now running ", spName, ", which is present on ", nS, " sites"))

        # add the data for the species of interest
        Cmodel$setData(spDat)

        # finish initialization
        spInits <- list(z = Z,
                        lam.0 = logit(dataSumm$stats$naiveOcc[sp] * .99)) # to avoid numeric problem

        # issue occurs here
        Cmodel$setInits(spInits)

        # test whether the model is fully initialised
        if(is.na(Cmodel$calculate())) {
          print(Cmodel$initializeInfo())
          stop("model not fully initialized")
          }

        # and now we can use $run on the compiled model object.
        samplesList <- list()
        if(mon2) samplesList2 <- list() else samplesList2 <- NULL

        for(i in 1:n.chain){
          CoccMCMC$run(niter = n.iter,
                       nburnin = n.burn,
                       chain = i,
                       thin = n.thin,
                       thin2 = 2 * n.thin, # for the annual parameters
                       reset = TRUE)
          samplesList[[i]] <- as.matrix(CoccMCMC$mvSamples)
          if(mon2) samplesList2[[i]] <- as.matrix(CoccMCMC$mvSamples2)

        }
        samplesList <- coda::as.mcmc.list(lapply(samplesList, as.mcmc))
        if(mon2) samplesList2 <- coda::as.mcmc.list(lapply(samplesList2, as.mcmc))

        return(list(fixed = samplesList, annual = samplesList2))
      }