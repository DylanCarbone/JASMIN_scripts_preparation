modified_summariseMSI <- 
function (label, plotType, indicator, method, plot = TRUE, writePlot = TRUE)  # removed year arguments, years chosen in calcMSI instead
{
  n <- max(indicator$Summary$Species_Number)
  nSpec <- indicator$Summary$Species_Number
  years <- indicator$Summary$year
  if (plotType == "indicator") {
    p1 <- ggplot(data = indicator$Summary[indicator$Summary$year %in% years,], aes(x = years, y = indicator)) + 
      geom_ribbon(data = indicator$Summary[indicator$Summary$year %in% years,], aes(ymax = upper, 
                                   ymin = lower), fill = "grey80") + 
      geom_line() + geom_point() + theme_linedraw() + ylab("Occupancy index") + 
      xlab("") + ylim(c(0, 180)) + ggtitle(label) + annotate("text", 
                                                             x = 1985, y = 30, label = paste(n, "species"))
    if (method == "lambda") {
      st <- indicator$st$species_assessment$category
      st <- data.frame(st, rep(as.factor(1), length(st)))
      colnames(st) <- c("val", "type")
      lt <- indicator$lt$species_assessment$category
      lt <- data.frame(lt, rep(as.factor(2), length(lt)))
      colnames(lt) <- c("val", "type")
      dat <- rbind(st, lt)
      p2 <- ggplot(dat, aes(x = factor(type), fill = forcats::fct_rev(val))) + 
        geom_bar(position = "fill") + theme_linedraw() + 
        ylab("Proportion of species") + xlab("") + scale_x_discrete(labels = c("Short term", 
                                                                               "Long term")) + guides(fill = guide_legend(title = ""))
      IndPlot <- grid.arrange(p1, p2, ncol = 2)
      print(IndPlot)
      if(writePlot){
        ggsave(IndPlot, filename = paste0("ind_", label, ".png"), height = 3, width = 9,
             units = "in", dpi = 300)
      }
    }
    else {
      print(p1)
    }
  }
  else if (plotType == "nSpecies") {
    p1 <- ggplot(data = NULL, aes(x = years, y = ind$summary$Species_Number)) + 
      geom_line() + geom_point() + theme_linedraw() + ylab("Number of species contributing") + 
      xlab("") + ggtitle(label)
    print(p1)
  }
  else {
    width <- summary$upper - summary$lower
    p1 <- ggplot(data = NULL, aes(y = width, x = ind$summary$Species_Number)) + 
      geom_point() + theme_linedraw() + ylab("Width of credible interval") + 
      xlab("Number of species contributing") + ggtitle(label)
    print(p1)
  }
  if (!plot) {
    indicator$Summary$year <- indicator$Summary$year + (minYear - 
                                                          1)
    out1 <- indicator$Summary
    changeLT <- data.frame(table(indicator$MetaData$species_change$category))
    rawLT <- indicator$MetaData$species_change[, 1]
    L <- length(rawLT)
    changeST <- data.frame(table(indicator$st$species_assessment$category))[, 
                                                                            2]
    rawST <- indicator$st$species_assessment[, 1]
    rawST <- c(rawST, rep(NA, (L - length(rawST))))
    changeFinal <- data.frame(table(indicator$final$species_assessment$category))[, 
                                                                                  2]
    rawFinal <- indicator$final$species_assessment[, 1]
    rawFinal <- c(rawFinal, rep(NA, (L - length(rawFinal))))
    change <- data.frame(changeLT, changeST, changeFinal)
    colnames(change) <- c("Trend", "Long term", "Short term", 
                          "Final year")
    raw <- data.frame(rawLT, rawST, rawFinal)
    colnames(raw) <- c("Long term", "Short term", "Final year")
    out <- list(out1, change, raw)
  }
  else {
    out <- p1
  }
}