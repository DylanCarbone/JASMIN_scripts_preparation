# index = 1;
# modPath = "/data-s3/occmods/"; 
# metaPath = "/data-s3/metadata/";
# ver = "most_recent";
# indicator = "all"; # original function doesn't work with "priority": there is no column in speciesInfo for UK or GB, so crashes in first line of applySamp. Kath has modified function to fix.
# region = "UK"; # doesn't work with "UK" but does with "GB":  Ants metadata doesn't have "UK" in region_aggs. Is this because no trends available for ants in NI? (there are none tagged in the region section of mr - only ENG, SCO, and WAL)
# nSamps = 999;
# minObs = 1;
# scaleObs = "global"; 
# write = TRUE; 
# outPath = "~/UKBI2021/filtered_D1c/";
# group = taxa;
# t0 = 1970; 
# tn = 2020;
# drop = FALSE
modified_createRoster <- 
  function (index, modPath = "/data-s3/occmods/", metaPath, ver = "most_recent", 
            group, indicator, region, nSamps = 999, minObs = 50, scaleObs = "global", 
            write, outPath, speciesToKeep = NULL, drop = TRUE, clipBy = "group", 
            t0, tn) 
  {
    if ("most_recent" %in% ver) {
      mr_files <- list.files("/data-s3/most_recent_meta")
      mr_files_ver <- gsub("metadata_", "", mr_files)
      mr_files_ver <- as.numeric(gsub(".csv", "", mr_files_ver))
      mr <- read.csv(paste0("/data-s3/most_recent_meta/", mr_files[which.max(mr_files_ver)]), 
                     stringsAsFactors = FALSE)
      tdf <- data.frame(ver_orig = ver, group = group) # LB changed taxa to group
      mr <- mr[mr$group %in% tdf$group & mr$most_recent == TRUE &  
                 mr$data_type == "occmod_outputs", ]
      mr_tdf <- merge(tdf, mr, by = "group") # LB changed taxa to group
      ver <- ifelse(mr_tdf$ver == "most_recent", mr_tdf$dataset_name, 
                    mr_tdf$ver)
    }
    if (all(region %in% c("GB", "UK", "ENGLAND", "SCOTLAND", 
                          "WALES", "NORTHERN_IRELAND")) == FALSE) {
      warning("Warning: not all regions match either GB, UK, ENGLAND, SCOTLAND, WALES, or NORTHERN_IRELAND is this what you expected?")
    }
    if(!is.null(speciesToKeep)){  #LB ifelse statement in the dataframe below creates a vector of logicals, not one
      speciesToKeep <- as.character(paste(speciesToKeep,  # LB
                                          collapse = ","))# LB
    }
    else speciesToKeep <- NA
    df <- data.frame(index = index, modPath = modPath, datPath = paste0("/data-s3/occmods/", 
                                                                        group, "/", ver, "/"), metaPath = metaPath, ver = ver, 
                     group = group, indicator = indicator, region = region, 
                     nSamps = as.numeric(nSamps), minObs = minObs, scaleObs = scaleObs, 
                     write = write, outPath = outPath, speciesToKeep = speciesToKeep, drop = drop, clipBy = clipBy, #LB
                     t0 = t0, tn = tn, stringsAsFactors = FALSE)
    roster <- split(df, seq(nrow(df)))
  }