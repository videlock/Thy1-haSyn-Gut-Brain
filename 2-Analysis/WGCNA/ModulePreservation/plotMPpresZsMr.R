plotMPpresZsMr<-function(mpdat,ref,test,filePrefix){
  modColors = rownames(mpdat$preservation$observed[[ref]][[test]])
  moduleSizes = mpdat$preservation$Z[[ref]][[test]][, 1];
  plotMods = !(modColors %in% c("grey", "gold"));
  text = modColors[plotMods];
  plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2],
                   mp$preservation$Z[[ref]][[test]][, 2])
  # Main titles for the plot
  mains = c("Preservation Median rank", "Preservation Zsummary");
  
  pdf(paste0(filePrefix,"_modulePreservation-Zsummary-medianRank.pdf"),
      wi=10, h=5)
  par(mfrow = c(1,2))
  par(mar = c(4.5,4.5,2.5,1)) 
  
  for (p in 1:2){
    min = min(plotData[, p], na.rm = TRUE)
    max = max(plotData[, p], na.rm = TRUE); # Adjust ploting ranges appropriately 
    if (p==2){
      if (min > -max/10) min = -max/10
      ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
    } else
      ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
    plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1,
         bg = modColors[plotMods], pch = 21,
         main = mains[p],
         cex = 2.4, ylab = mains[p], xlab = "Module size", log = "x",
         ylim = ylim, xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
    
    # For Zsummary, add threshold lines
    if (p==2) {
      labelPoints(x = moduleSizes[plotMods], y = plotData[plotMods, p], 
                  text, cex = .5, offs = 0.08);
      abline(h=0)
      abline(h=2, col = "blue", lty = 2)
      abline(h=10, col = "darkgreen", lty = 2)
    } 
  }
  dev.off()
}
