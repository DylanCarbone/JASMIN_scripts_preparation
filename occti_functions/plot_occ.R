# Plot index
#==================================================================


plot_occ <- function(output){

z$psiA_Ub <- sapply(z$psiA_U,function(x){min(x,1)})
z$psiA_Lb <- sapply(z$psiA_L,function(x){max(x,0)})

  jpeg(filename=paste("./Output/Plots/",spp, "_occ_index_",minyear,"_",maxyear,plottitle,res,".jpg", sep=""), quality=100, width=800)
  print(ggplot(z, aes_string(x = "Year", y = "psiA")) +
          theme_bw() +
          geom_ribbon(aes(ymin = psiA_Lb,
                          ymax = psiA_Ub), alpha = 0.2) +
          geom_line(size = 1, col = "black") +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),text=element_text(size=13))+
          ylab("Occupancy index") +
          xlab("Year"))
  dev.off()

  jpeg(filename=paste("./Output/Plots/",spp, "_occ_index_",minyear,"_",maxyear,"_01",plottitle,res,".jpg", sep=""), quality=100, width=800)
  print(ggplot(z, aes_string(x = "Year", y = "psiA")) +
          theme_bw() +
          geom_ribbon(aes(ymin = psiA_Lb,
                          ymax = psiA_Ub), alpha = 0.2) +
          geom_line(size = 1, col = "black") +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),text=element_text(size=13))+
          ylab("Occupancy index") +
          xlab("Year") + ylim(0,1))
  dev.off()

  if(boot){
    jpeg(filename=paste("./Output/Plots/",spp, "_occ_index_",minyear,"_",maxyear,"_unbounded_boot",plottitle,res,".jpg", sep=""), quality=100, width=800)
    print(ggplot(z, aes_string(x = "Year", y = "psiA")) +
            theme_bw() +
            geom_ribbon(aes(ymin = psiA_L_boot,
                            ymax = psiA_U_boot), alpha = 0.2) +
            geom_line(size = 1, col = "black") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),text=element_text(size=13))+
            ylab("Occupancy index") +
            xlab("Year"))
    dev.off()

    jpeg(filename=paste("./Output/Plots/",spp, "_occ_index_",minyear,"_",maxyear,"_boot",plottitle,res,".jpg", sep=""), quality=100, width=800)
    print(ggplot(z, aes_string(x = "Year", y = "psiA")) +
            theme_bw() +
            geom_ribbon(aes(ymin = psiA_L_boot,
                            ymax = psiA_U_boot), alpha = 0.2) +
            geom_line(size = 1, col = "black") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),text=element_text(size=13))+
            ylab("Occupancy index") +
            xlab("Year") + ylim(0,1))
    dev.off()
  }
  if(plotrecords){
    recs <- obdata_sp[, .N, by=Year]

    jpeg(filename=paste("./Output/Plots/",spp, "_nrecords_",minyear,"_",maxyear,plottitle,res,".jpg", sep=""), quality=100, width=800)
    print(ggplot(recs, aes_string(x = "Year", y = "N")) +
            theme_bw() +
            geom_line(size = 1, col = "black") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),text=element_text(size=13))+
            ylab("Number of records") +
            xlab("Year"))
    dev.off()

    y20 <- min(subset(recs, N >= 20)$Year)
    y30 <-  min(subset(recs, N >= 30)$Year)

    jpeg(filename=paste("./Output/Plots/",spp, "_occ_index_",minyear,"_",maxyear,"_01_recs",plottitle,res,".jpg", sep=""), quality=100, width=800)
    print(ggplot(z, aes_string(x = "Year", y = "psiA")) +
            theme_bw() +
            geom_ribbon(aes(ymin = psiA_Lb,
                            ymax = psiA_Ub), alpha = 0.2) +
            geom_line(size = 1, col = "black") +
            geom_vline(xintercept=y20, linetype="dashed",color="blue", size=1)+
            geom_vline(xintercept=y30, linetype="dashed",color="green", size=1)+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),text=element_text(size=13))+
            ylab("Occupancy index") +
            xlab("Year") + ylim(0,1))
    dev.off()


}}


