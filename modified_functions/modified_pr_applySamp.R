modified_pr_applySamp <-
  function (roster, parallel = TRUE, sample = TRUE){
    data(speciesInfo)
    if (roster$indicator %in% c("priority", "pollinators")) { #LB
      if(roster$indicator == "priority"){
      
      ## KATTUR replaced this ...
      
      # keepInds <- which(!is.na(speciesInfo[, roster$region]))
      # keep <- c(as.character(speciesInfo$Species[keepInds]), 
      #           as.character(speciesInfo$concept[keepInds]))
      # keep <- keep[-which(is.na(keep))]
      
      
      ## With this ...
      
      if (!roster$region %in% c("WAL", "WALES", "SCO", "SCOTLAND",
                                "ENG", "ENGLAND", "NIR","NORTHERN_IRELAND", "GB", "UK")){ #LB added "GB" and "UK"
        
        stop("Priority species region not recognised")
      }
      
      if(roster$region == "WAL" | roster$region == "WALES") {pr_region <- "WALES"}
      if(roster$region == "ENG" | roster$region == "ENGLAND") {pr_region <- "ENGLAND"}
      if(roster$region == "SCO" | roster$region == "SCOTLAND") {pr_region <- "SCOTLAND"}
      if(roster$region == "NIR" | roster$region == "NORTHERN_IRELAND") {pr_region <- "NORTHERN.IRELAND"}
      if(roster$region %in% c("UK")) {pr_region <- c("WALES","ENGLAND", "SCOTLAND","NORTHERN.IRELAND" )} #LB
      if(roster$region %in% c("GB")) {pr_region <- c("WALES","ENGLAND", "SCOTLAND")}
      
      
      speciesInfo_group <- speciesInfo[speciesInfo$Group %in% roster$group,]
      
      if (!any(!is.na(speciesInfo_group[, pr_region,]))) {
        stop("No priority species in this group")
      } else {
        #LB replaced
        #keepInds <- which(!is.na(speciesInfo_group[, pr_region,]))
        # with
        keepInds <- sapply(1:nrow(speciesInfo_group), function(i) any(speciesInfo_group[i,pr_region] == "Y"))  
        # as before from here
        keep <- c(as.character(speciesInfo_group$Species[keepInds]), 
                  as.character(speciesInfo_group$concept[keepInds]))
        keep <- unique(keep)
        }
      }
      if(roster$indicator == "pollinators"){
        keep <- wrappeR:::sampSubset("pollinators", inPath = roster$metaPath)
        }
      modFilePath <- file.path(roster$modPath, roster$group, 
                                 "occmod_outputs", roster$ver)
      modFiles_rdata <- list.files(modFilePath, pattern = ".rdata")
      modFiles_rds <- list.files(modFilePath, pattern = ".rds")
      if (length(modFiles_rdata) == 0 & length(modFiles_rds) ==  0) {
            stop("Model files must be either .rds or .rdata")
          } else {
           if (length(modFiles_rdata) == 0) {
             filetype <- "rds"
             modFiles <- modFiles_rds
           } else {
             filetype <- "rdata"
             modFiles <- modFiles_rdata 
             }
         }
      keep_iter <- gsub(pattern = paste0("\\.", filetype), 
                           repl = "", modFiles) 
      first_spp <- keep_iter[[1]]
      if (substr(first_spp, (nchar(first_spp) + 1) - 2,
                    nchar(first_spp)) %in% c("_1", "_2", "_3")) {
      keep_name <- keep[tolower(keep) %in% tolower(stringr::str_replace(keep_iter,"(?<![A-z][A_z][A-z])_.*", ""))]
      keep_concept<- keep[tolower(keep) %in% tolower(stringr::str_replace(keep_iter,"(?<![a-z])_.*", ""))]
      keep <- c(keep_concept, keep_name)
      keep <- unique(keep)
      keep <- keep[!is.na(keep)]
      keep_iter <- keep_iter[(tolower(stringr::str_replace(keep_iter,"(?<![a-z])_.*", ""))
                                    %in% tolower(keep)) |
                                       (tolower(stringr::str_replace(keep_iter,"(?<![A-z][A_z][A-z])_.*", ""))
                                      %in% tolower(keep))]
          }  else {
           keep <- keep_iter[tolower(keep_iter) %in% tolower(keep)]
           keep_iter <- NULL
           keep <- unique(keep)
          }
      }
    
    if (roster$indicator == "all"){ ## Added by KATTUR, modified by LB
      modFilePath <- file.path(roster$modPath, roster$group, 
                               "occmod_outputs", roster$ver)
      modFiles_rdata <- list.files(modFilePath, pattern = ".rdata")
      modFiles_rds <- list.files(modFilePath, pattern = ".rds")
      if (length(modFiles_rdata) == 0 & length(modFiles_rds) == 
          0) 
        stop("Model files must be either .rds or .rdata")
      else {
        if (length(modFiles_rdata) == 0) {
          filetype <- "rds"
          modFiles <- modFiles_rds
        }
        else {
          filetype <- "rdata"
          modFiles <- modFiles_rdata
        }
      }
      keep_iter <- gsub(pattern = paste0("\\.", filetype), 
                        repl = "", modFiles)
      first_spp <- keep_iter[[1]]
      if (substr(first_spp, (nchar(first_spp) + 1) - 2, nchar(first_spp)) %in% 
          c("_1", "_2", "_3")) {
        keep <- gsub("(.*)_\\w+", "\\1", keep_iter)
        keep <- gsub("(.*)_\\w+", "\\1", keep)
        keep <- unique(keep)
      }
      else {
        keep <- keep_iter
        keep_iter <- NULL
      }
    } ## Added by KATTUR
    
    
    if (!is.na(roster$speciesToKeep)) {
      speciesToKeep <- unique(unlist(strsplit(roster$speciesToKeep, 
                                       ",")))
      
      # notFound <- speciesToKeep[!tolower(speciesToKeep) %in% #LB
      #                             tolower(keep)]             #LB
      notFound <- setdiff(speciesToKeep, speciesInfo$Species)  #LB
      if (length(notFound) > 0) {
        warning(paste("some species on your \"speciesToKeep\" list were not found in the data:", 
                      paste(notFound, collapse = ", ")))
      }
      # lookup concept codes, in case this is how some models are named (Bryophytes are a mix)
      spCon <- speciesInfo[speciesInfo$Species %in% speciesToKeep & 
                             speciesInfo$Group %in% roster$group,c("concept", "Species")]
      keep <- intersect(keep, c(as.character(spCon$concept), as.character(spCon$Species)))
  }
    
    
    
    if (roster$drop == TRUE) {
      drop <- which(!is.na(speciesInfo$Reason_not_included) & 
                      speciesInfo$Reason_not_included != "Didn't meet criteria")
      drop <- c(as.character(speciesInfo$Species[drop]), as.character(speciesInfo$concept[drop]))
      keep <- keep[!keep %in% drop]
    }
    
    
    
    if (sample == TRUE) 
      out <- modified_tempSampPost(indata = paste0(roster$modPath, roster$group, 
                                                   "/occmod_outputs/", roster$ver, "/"), keep = keep, 
                                   keep_iter = keep_iter, output_path = NULL, region = roster$region, 
                                   sample_n = roster$nSamps, group_name = roster$group, 
                                   combined_output = TRUE, write = FALSE, minObs = roster$minObs, 
                                   scaleObs = roster$scaleObs, t0 = roster$t0, tn = roster$tn, 
                                   parallel = parallel, filetype = filetype)
    
    
    else out <- getA(indata = paste0(roster$modPath, roster$group, 
                                     "/occmod_outputs/", roster$ver, "/"), keep = keep, REGION_IN_Q = paste0("a_", 
                                                                                                             roster$region), group_name = roster$group, combined_output = TRUE, 
                     write = FALSE, minObs = roster$minObs, t0 = roster$t0, 
                     tn = roster$tn, parallel = parallel)
    samp_post <- out[[1]]
    samp_post$species <- tolower(samp_post$species)
    meta <- out[[2]]
    meta[, 1] <- tolower(meta[, 1])
    if (roster$write == TRUE) {
      comb <- list(samp_post, meta)
      save(comb, file = paste0(roster$outPath, roster$group, 
                               "_", roster$indicator, "_", roster$region, "_samp.rdata"))
    }
    return(list(samp_post = samp_post, meta = meta, indicator = roster$indicator, 
                group_name = roster$group, region = roster$region, clipBy = roster$clipBy, 
                minObs = roster$minObs))
  }