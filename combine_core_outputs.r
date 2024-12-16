    names(yearEff) <- dimnames(transferable_outputs$dataSumm$occMatrix)[[1]][species_indices][core_i]
    attr(yearEff, "modelCode") <- transferable_outputs$model$getCode()
    saveRDS(yearEff, file = paste0("yearEff_node_", node_i, ".rds"))