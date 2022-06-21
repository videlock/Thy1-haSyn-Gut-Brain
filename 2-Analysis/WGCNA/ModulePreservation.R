
options(stringsAsFactors = FALSE);
library(WGCNA);
allowWGCNAThreads() 
library(flashClust);
library(gplots)
library(RColorBrewer)

# plot function---------
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

# colon specific/consensus----------

# create input----------

load("../cons/inputData.rda")
datExprCons<-multiExpr$dc$data
load("../cons/modules.rda")
colorsCons<-net$colors

nSets = 2
setLabels = c("ColonSpecific", "ColonConsensus");
multiExpr = list()
colorList = list()

stdGenes = NULL

load("../dc/processed_data/inputData.rda")
genes = colnames(datExpr)
stdGenes = genes
multiExpr$ColonSpecific$data<-datExpr
load("../dc/modules.rda")
colorList$ColonSpecific<-moduleColors

genes= colnames(datExprCons)
stdGenes = intersect(genes, stdGenes)

multiExpr$ColonConsensus$data<-datExprCons
colorList$ColonConsensus<-colorsCons

# Restrict data to common genes
for(set in setLabels){
  keep=colnames(multiExpr[[set]]$data)%in%stdGenes
  multiExpr[[set]]$data<-multiExpr[[set]]$data[,keep]
  colorList[[set]]<-colorList[[set]][keep]
}

lapply(multiExpr, lapply, dim)


save(multiExpr,colorList,file = "input_mp_dcSpvC.rda")

# run pres-------
rm(list = ls())
source("helperFunctions.R")

load("input_mp_dcSpvC.rda")
nSets<-length(multiExpr)
setLabels=names(multiExpr)

mp = modulePreservation(multiExpr, colorList,networkType = "signed",
                        corFnc = "bicor",
                        corOptions = c("maxPOutliers = 0.05",
                                       "maxPOutliers = 0.05"),
                        referenceNetworks = c(1,2),
                        loadPermutedStatistics = FALSE,
                        savePermutedStatistics = TRUE,
                        permutedStatisticsFile = "permutedStats_dcSpvC.rda",
                        nPermutations = 100,
                        verbose = 3)

save(mp,colorList,multiExpr,nSets,setLabels,file = "mpdata_dcSpvC.rda")


# calculate statistics ----------------
load("mpdata_dcSpvC.rda")
source("helperFunctions.R")

library(clusterRepro)
library(impute)

eigengenes = list();
for (set in 1:nSets)
{
  eigengenes[[set]] = multiSetMEs(multiExpr, universalColors = colorList[[set]], 
                                  excludeGrey = TRUE);
  for (ss in 1:nSets)
  {
    rownames(eigengenes[[set]][[ss]]$data) = rownames(multiExpr[[ss]]$data);
  }
}

ref = 1 
test = 2 
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], 
                 mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], 
               mp$preservation$Z[[ref]][[test]][, -1]);
print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2));
# Compare preservation to quality:
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

write.csv(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2),
          "mpstats_dcSpvC_Spref.csv")


plotMPpresZsMr(mpdat = mp,ref = 1,test = 2,filePrefix = "dcSpPresInCons")

ref = 2
test = 1
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], 
                 mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], 
               mp$preservation$Z[[ref]][[test]][, -1]);
print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2));
# Compare preservation to quality:
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

write.csv(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2),
          "mpstats_dcSpvC_Cpref.csv")


plotMPpresZsMr(mpdat = mp,ref = 2,test = 1,filePrefix = "dcConsPresInSp")


# summary table--------

# This variable will contain the summary table
summaryTable = NULL
# Loop over all combinations of reference and tests sets
for (ref in 1:nSets) for (test in 1:nSets) if (ref!=test){
  modules = rownames(mp$preservation$Z[[ref]][[test]])
  nMods = length(modules)
  sizes = mp$preservation$Z[[ref]][[test]][, 1]
  acc = matrix(NA, nMods, 3)
  # if(test!=2){
  acc[match(rownames(mp$accuracy$observed[[ref]][[test]]), modules), ] = mp$accuracy$observed[[ref]][[test]][, -1, drop = FALSE]
  colnames(acc) = colnames(mp$accuracy$observed[[ref]][[test]])[-1]
  accZ = mp$accuracy$Z[[ref]][[test]][, -1, drop = FALSE]
  acc.log.p = mp$accuracy$log.p[[ref]][[test]][, -1, drop = FALSE]
  acc.log.pBonf = mp$accuracy$log.pBonf[[ref]][[test]][, -1, drop = FALSE]
  # }else{
  accZ = matrix(NA, nMods, 3)
  acc.log.p = matrix(NA, nMods, 3)
  acc.log.pBonf = matrix(NA, nMods, 3)
  colnames(acc) = colnames(mp$accuracy$observed[[1]][[2]])[-1]
  colnames(accZ) = colnames(mp$accuracy$Z[[1]][[2]])[-1]
  colnames(acc.log.p) = colnames(mp$accuracy$log.p[[1]][[2]])[-1]
  colnames(acc.log.pBonf) = colnames(mp$accuracy$log.pBonf[[1]][[2]])[-1]
  # }
  ##Table of results for this reference-test combination
  tab = cbind(referenceSet = rep(setLabels[ref], nMods),
              testSet = rep(setLabels[test], nMods),
              moduleLabel = modules,
              moduleSize = sizes,
              mp$quality$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$observed[[ref]][[test]][, -1, drop = FALSE],
              acc,
              mp$referenceSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              accZ,
              acc.log.p,
              acc.log.pBonf,
              mp$referenceSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$referenceSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$referenceSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE]
  )
  # Add the table to the main table.
  if (is.null(summaryTable)){
    summaryTable = tab 
  }else{
    summaryTable = rbind(summaryTable, tab)
  }
}
rm(acc,acc.log.p,acc.log.pBonf,accZ,tab)
# Save the table in csv format.
write.table(summaryTable, file = "completeResults_dcSpvC.csv", 
            row.names = FALSE, sep = ",", quote = FALSE)


# scaled expression------
adjList=list()
adjList$ColonSpecific<-adjacency(multiExpr$ColonSpecific$data, power = 8, type = "signed",
                                 corFnc="bicor",corOptions = list(maxPOutliers =0.05))

adjList$ColonConsensus<-adjacency(multiExpr$ColonConsensus$data, power = 8, type = "signed",
                                  corFnc="bicor",corOptions = list(maxPOutliers =0.05))

IMconlist=list()
for(set in setLabels){
  IMconlist[[set]]<-apply(adjList[[set]],1,sum)
}

scaledIMconlist=list()
for(set in setLabels){
  scaledIMconlist[[set]]<-IMconlist[[set]]/max(IMconlist[[set]])
}

rm(IMconlist)

scaledExprList=list()
for(set in setLabels){
  meanExpr<-apply(multiExpr[[set]][["data"]],2,mean)
  scaledExprList[[set]]<-meanExpr/max(meanExpr)
}

lapply(scaledIMconlist,median)



minconnections=.3
goodIMC<-lapply(scaledIMconlist,function(x) return(ifelse(x>minconnections,1,0)))

sum(goodIMC$ColonSpecific) #all 9052
sum(goodIMC$ColonConsensus) #all 9052

keep=goodIMC$ColonSpecific>0&goodIMC$ColonConsensus>0
table(keep)


adjList.rest<-lapply(adjList,function(x) return(x[keep,keep]))
rm(adjList)

TOMdistList<-lapply(adjList.rest,TOMdist)

dendroList<-lapply(TOMdistList, function(x) flashClust(as.dist(x), method = "a"))

meList<-list()
for(refset in setLabels){
  for(set in setLabels){
    meList[[refset]] = moduleEigengenes(multiExpr[[set]]$data[,keep], 
                                        colorList[[refset]][keep])$eigengenes
  }
}

# heatmap---------
setPairs=NULL
for (set2 in setLabels){
  for (set1 in setLabels){
    if (set1!=set2){
      setpair=paste(set1,set2,sep = "_")
      if (!paste(set2,set1,sep ="_")%in%setPairs[,1]){
        setPairs<-rbind(setPairs,c(setpair,set1,set2))
      }
    }
  }
}
setPairs<-as.data.frame(setPairs)
names(setPairs)<-c("setpair","set1","set2")
setNames<-list(ColonConsensus="Colon-Striatum Consensus",
               ColonSpecific="Colon Specific")
setModFiles<-list(ColonConsensus="../cons/modtab.rds",
                  ColonSpecific="../dc/modtab.rds")
overlapList<-list()

for(i in 1:nrow(setPairs)){
  setpair=setPairs$setpair[i]
  set1=setPairs$set1[i]
  set2=setPairs$set2[i]
  setname1=setNames[[set1]]
  setname2=setNames[[set2]]
  modtab1<-readRDS(setModFiles[[set1]])
  modtab2<-readRDS(setModFiles[[set2]])
  overlap = overlapTable(colorList[[set1]][keep],
                         colorList[[set2]][keep],
                         levels1 =modtab1$moduleColor,
                         levels2 = modtab2$moduleColor)
  numMat = -log10(overlap$pTable);
  numMat[numMat >50] = 50;
  # Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of # counts and corresponding p-values.
  textMat = paste0(overlap$countTable, "\n", signif(overlap$pTable, 2));
  dim(textMat) = dim(numMat)
  # Additional information for the plot. These will be used shortly.
  # xLabels = paste("M", sort(unique(colorList[[set2]])));
  # yLabels = paste("M", sort(unique(colorList[[set1]])));
  xLabels=rownames(modtab2)
  yLabels=rownames(modtab1)
  
  # xSymbols = paste(sort(unique(colorList[[set2]])), ": ", 
  #                  table(colorList[[set2]][keep]), sep = "")
  # ySymbols = paste(sort(unique(colorList[[set1]])), ": ", 
  #                  table(colorList[[set1]][keep]), sep = "")
  
  xSymbols=paste0(modtab2$moduleColor,"(",modtab2$Size,")")
  ySymbols=paste0(modtab1$moduleColor,"(",modtab1$Size,")")
  
  xSymbols2=paste0(modtab2$ShortLabel,"(",modtab2$Size,")")
  ySymbols2=paste0(modtab1$ShortLabel,"(",modtab1$Size,")")
  
  
  
  
  overlapList[[setpair]]<-list(set1=set1,
                               set2=set2,
                               setname1=setname1,
                               setname2=setname2,
                               overlap=overlap,
                               numMat=numMat,
                               textMat=textMat,
                               xLabels=xLabels,
                               yLabels=yLabels,
                               xSymbols=xSymbols,
                               ySymbols=ySymbols,
                               xSymbols2=xSymbols2,
                               ySymbols2=ySymbols2
  )
  
}

setpair=setPairs[1,1]
pdat=overlapList[[setpair]]
filename=paste0(setpair,"_dendrosAndTable.pdf")
cex.main = 1.2
cex.textmatrix=.6
cex.matrixlab=.9

pdf(filename,w = 12, h = 6)
fp = TRUE; layout(matrix(c(1,2,5, 3,4,5), 3, 2), 
                  heights = c(1.5, 0.5, 8), widths = c(1, 1)); par(mgp = c(3, 1, 0))
plotDendroAndColors(dendroList[[pdat$set1]], 
                    cbind(colorList[[pdat$set1]][keep], 
                          colorList[[pdat$set2]][keep]),
                    c("Specific","Consensus"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2), addGuide = FALSE,
                    main = paste0("A. ",pdat$setname1," gene dendrogram and module colors"),
                    cex.main = cex.main, dendroLabels = FALSE, hang = 0.03,
                    autoColorHeight = FALSE, colorHeight = 0.5,cex.axis=.6,
                    cex.colorLabels = 0.7, abHeight = 0.95);par(mgp = c(3, 1, 0));
plotDendroAndColors(dendroList[[pdat$set2]], 
                    cbind(colorList[[pdat$set1]][keep], 
                          colorList[[pdat$set2]][keep]),
                    c("Specific","Consensus"),
                    # c(paste(pdat$set1,"modules"),paste(pdat$set2,"modules")), 
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2), addGuide = FALSE,
                    main = paste0("B. ",pdat$setname2," gene dendrogram and module colors"),
                    cex.main = cex.main, dendroLabels = FALSE, hang = 0.03,
                    autoColorHeight = FALSE, colorHeight = 0.5,cex.axis=.6,
                    cex.colorLabels = 0.7, abHeight = 0.95);
par(mar = c(8, 7.5, 1.5, .7));
labeledHeatmap(Matrix = pdat$numMat,
               xLabels = pdat$xLabels,xSymbols = pdat$xSymbols,
               yLabels = pdat$yLabels,ySymbols = pdat$ySymbols,
               colorLabels = TRUE, colors = blueWhiteRed(200)[100:200],
               textMatrix = pdat$textMat, cex.text = cex.textmatrix,
               setStdMargins = FALSE, cex.lab = cex.matrixlab,
               xColorWidth = 0.01,xLabelsAngle = 90,yColorWidth = 0.01,
               verticalSeparator.col = "lightgrey",
               verticalSeparator.interval = 1,
               horizontalSeparator.interval = 1,
               horizontalSeparator.col = "lightgrey",
               main = paste0("C. ",pdat$setname1," modules (rows) vs. ", pdat$setname2, " modules (columns)"),
               cex.main = cex.main);dev.off()

save(colorList,dendroList,
     modtab1,modtab2,
     meList,multiExpr,overlapList,scaledExprList,
     scaledIMconlist,setPairs,goodIMC,statsObs,summaryTable,cex.main,
     cex.matrixlab,cex.textmatrix,
     keep,setLabels,file = "mp_overlap_dcSpvC.rda")



# striatum specific/consensus ---------------

# create input----------

load("../cons/inputData.rda")
datExprCons<-multiExpr$str$data
load("../cons/modules.rda")
colorsCons<-net$colors

nSets = 2
setLabels = c("StriatumSpecific", "StriatumConsensus");
multiExpr = list()
colorList = list()

stdGenes = NULL

load("../str/processed_data/inputData.rda")
genes = colnames(datExpr)
stdGenes = genes
multiExpr$StriatumSpecific$data<-datExpr
load("../str/modules.rda")
colorList$StriatumSpecific<-moduleColors

genes= colnames(datExprCons)
stdGenes = intersect(genes, stdGenes)

multiExpr$StriatumConsensus$data<-datExprCons
colorList$StriatumConsensus<-colorsCons

# Restrict data to common genes
for(set in setLabels){
  keep=colnames(multiExpr[[set]]$data)%in%stdGenes
  multiExpr[[set]]$data<-multiExpr[[set]]$data[,keep]
  colorList[[set]]<-colorList[[set]][keep]
}

lapply(multiExpr, lapply, dim)


save(multiExpr,colorList,file = "input_mp_strSpvC.rda")

# run pres-------
rm(list = ls())

load("input_mp_strSpvC.rda")
nSets<-length(multiExpr)
setLabels=names(multiExpr)

mp = modulePreservation(multiExpr, colorList,networkType = "signed",
                        corFnc = "bicor",
                        corOptions = c("maxPOutliers = 0.01",
                                       "maxPOutliers = 0.01"),
                        referenceNetworks = c(1,2),
                        loadPermutedStatistics = FALSE,
                        savePermutedStatistics = TRUE,
                        permutedStatisticsFile = "permutedStats_strSpvC.rda",
                        nPermutations = 100,
                        verbose = 3)

save(mp,colorList,multiExpr,nSets,setLabels,file = "mpdata_strSpvC.rda")


# calculate statistics ----------------
load("mpdata_strSpvC.rda")
source("helperFunctions.R")

library(clusterRepro)
library(impute)

eigengenes = list();
for (set in 1:nSets)
{
  eigengenes[[set]] = multiSetMEs(multiExpr, universalColors = colorList[[set]], 
                                  excludeGrey = TRUE);
  for (ss in 1:nSets)
  {
    rownames(eigengenes[[set]][[ss]]$data) = rownames(multiExpr[[ss]]$data);
  }
}

ref = 1 
test = 2 
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], 
                 mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], 
               mp$preservation$Z[[ref]][[test]][, -1]);
print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2));
# Compare preservation to quality:
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

write.csv(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2),
          "mpstats_strSpvC_Spref.csv")


plotMPpresZsMr(mpdat = mp,ref = 1,test = 2,filePrefix = "strSpPresInCons")

ref = 2
test = 1
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], 
                 mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], 
               mp$preservation$Z[[ref]][[test]][, -1]);
print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2));
# Compare preservation to quality:
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

write.csv(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2),
          "mpstats_strSpvC_Cpref.csv")


plotMPpresZsMr(mpdat = mp,ref = 2,test = 1,filePrefix = "strConsPresInSp")


# summary table--------

# This variable will contain the summary table
summaryTable = NULL
# Loop over all combinations of reference and tests sets
for (ref in 1:nSets) for (test in 1:nSets) if (ref!=test){
  modules = rownames(mp$preservation$Z[[ref]][[test]])
  nMods = length(modules)
  sizes = mp$preservation$Z[[ref]][[test]][, 1]
  acc = matrix(NA, nMods, 3)
  # if(test!=2){
  acc[match(rownames(mp$accuracy$observed[[ref]][[test]]), modules), ] = mp$accuracy$observed[[ref]][[test]][, -1, drop = FALSE]
  colnames(acc) = colnames(mp$accuracy$observed[[ref]][[test]])[-1]
  accZ = mp$accuracy$Z[[ref]][[test]][, -1, drop = FALSE]
  acc.log.p = mp$accuracy$log.p[[ref]][[test]][, -1, drop = FALSE]
  acc.log.pBonf = mp$accuracy$log.pBonf[[ref]][[test]][, -1, drop = FALSE]
  # }else{
  accZ = matrix(NA, nMods, 3)
  acc.log.p = matrix(NA, nMods, 3)
  acc.log.pBonf = matrix(NA, nMods, 3)
  colnames(acc) = colnames(mp$accuracy$observed[[1]][[2]])[-1]
  colnames(accZ) = colnames(mp$accuracy$Z[[1]][[2]])[-1]
  colnames(acc.log.p) = colnames(mp$accuracy$log.p[[1]][[2]])[-1]
  colnames(acc.log.pBonf) = colnames(mp$accuracy$log.pBonf[[1]][[2]])[-1]
  # }
  ##Table of results for this reference-test combination
  tab = cbind(referenceSet = rep(setLabels[ref], nMods),
              testSet = rep(setLabels[test], nMods),
              moduleLabel = modules,
              moduleSize = sizes,
              mp$quality$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$observed[[ref]][[test]][, -1, drop = FALSE],
              acc,
              mp$referenceSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              accZ,
              acc.log.p,
              acc.log.pBonf,
              mp$referenceSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$referenceSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$referenceSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE]
  )
  # Add the table to the main table.
  if (is.null(summaryTable)){
    summaryTable = tab 
  }else{
    summaryTable = rbind(summaryTable, tab)
  }
}
rm(acc,acc.log.p,acc.log.pBonf,accZ,tab)
# Save the table in csv format.
write.table(summaryTable, file = "completeResults_strSpvC.csv", 
            row.names = FALSE, sep = ",", quote = FALSE)


# scaled expression------
adjList=list()
adjList$StriatumSpecific<-adjacency(multiExpr$StriatumSpecific$data, power = 12, type = "signed",
                                    corFnc="bicor",corOptions = list(maxPOutliers =0.01))

adjList$StriatumConsensus<-adjacency(multiExpr$StriatumConsensus$data, power = 12, type = "signed",
                                     corFnc="bicor",corOptions = list(maxPOutliers =0.01))

IMconlist=list()
for(set in setLabels){
  IMconlist[[set]]<-apply(adjList[[set]],1,sum)
}

scaledIMconlist=list()
for(set in setLabels){
  scaledIMconlist[[set]]<-IMconlist[[set]]/max(IMconlist[[set]])
}

rm(IMconlist)


lapply(scaledIMconlist,median)

minconnections=.3
goodIMC<-lapply(scaledIMconlist,function(x) return(ifelse(x>minconnections,1,0)))

sum(goodIMC$StriatumSpecific) 
sum(goodIMC$StriatumConsensus) 
# 7345

keep=goodIMC$StriatumSpecific>0&goodIMC$StriatumConsensus>0
table(keep)
# FALSE  TRUE 
# 1707  7345 

adjList.rest<-lapply(adjList,function(x) return(x[keep,keep]))
rm(adjList)

TOMdistList<-lapply(adjList.rest,TOMdist)

dendroList<-lapply(TOMdistList, function(x) flashClust(as.dist(x), method = "a"))

meList<-list()
for(refset in setLabels){
  for(set in setLabels){
    meList[[refset]] = moduleEigengenes(multiExpr[[set]]$data[,keep], 
                                        colorList[[refset]][keep])$eigengenes
  }
}

setPairs=NULL
for (set2 in setLabels){
  for (set1 in setLabels){
    if (set1!=set2){
      setpair=paste(set1,set2,sep = "_")
      if (!paste(set2,set1,sep ="_")%in%setPairs[,1]){
        setPairs<-rbind(setPairs,c(setpair,set1,set2))
      }
    }
  }
}
setPairs<-as.data.frame(setPairs)
names(setPairs)<-c("setpair","set1","set2")
setNames<-list(StriatumConsensus="Colon-Striatum Consensus",
               StriatumSpecific="Striatum Specific")



setModFiles<-list(StriatumConsensus="../cons/modtab.rds",
                  StriatumSpecific="../str/modtab-g.rds")
overlapList<-list()

for(i in 1:nrow(setPairs)){
  setpair=setPairs$setpair[i]
  set1=setPairs$set1[i]
  set2=setPairs$set2[i]
  setname1=setNames[[set1]]
  setname2=setNames[[set2]]
  modtab1<-readRDS(setModFiles[[set1]])
  modtab2<-readRDS(setModFiles[[set2]])
  overlap = overlapTable(colorList[[set1]][keep],
                         colorList[[set2]][keep],
                         levels1 =modtab1$moduleColor,
                         levels2 = modtab2$moduleColor)
  numMat = -log10(overlap$pTable);
  numMat[numMat >50] = 50;
  # Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of # counts and corresponding p-values.
  textMat = paste0(overlap$countTable, "\n", signif(overlap$pTable, 2));
  dim(textMat) = dim(numMat)
  xLabels=rownames(modtab2)
  yLabels=rownames(modtab1)
  
  xSymbols=paste0(modtab2$moduleColor,"(",modtab2$Size,")")
  ySymbols=paste0(modtab1$moduleColor,"(",modtab1$Size,")")
  
  
  
  
  overlapList[[setpair]]<-list(set1=set1,
                               set2=set2,
                               setname1=setname1,
                               setname2=setname2,
                               overlap=overlap,
                               numMat=numMat,
                               textMat=textMat,
                               xLabels=xLabels,
                               yLabels=yLabels,
                               xSymbols=xSymbols,
                               ySymbols=ySymbols
  )
  
}

setpair=setPairs[1,1]
pdat=overlapList[[setpair]]
filename=paste0(setpair,"_dendrosAndTable.pdf")
cex.main = 1.2
cex.textmatrix=.6
cex.matrixlab=.9

pdf(filename,w = 12, h = 6)
fp = TRUE; layout(matrix(c(1,2,5, 3,4,5), 3, 2), 
                  heights = c(1.5, 0.5, 8), widths = c(1, 1)); par(mgp = c(3, 1, 0))
plotDendroAndColors(dendroList[[pdat$set1]], 
                    cbind(colorList[[pdat$set1]][keep], 
                          colorList[[pdat$set2]][keep]),
                    c("Specific","Consensus"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2), addGuide = FALSE,
                    main = paste0("A. ",pdat$setname1," gene dendrogram and module colors"),
                    cex.main = cex.main, dendroLabels = FALSE, hang = 0.03,
                    autoColorHeight = FALSE, colorHeight = 0.5,cex.axis=.6,
                    cex.colorLabels = 0.7, abHeight = 0.95);par(mgp = c(3, 1, 0));
plotDendroAndColors(dendroList[[pdat$set2]], 
                    cbind(colorList[[pdat$set1]][keep], 
                          colorList[[pdat$set2]][keep]),
                    c("Specific","Consensus"),
                    # c(paste(pdat$set1,"modules"),paste(pdat$set2,"modules")), 
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2), addGuide = FALSE,
                    main = paste0("B. ",pdat$setname2," gene dendrogram and module colors"),
                    cex.main = cex.main, dendroLabels = FALSE, hang = 0.03,
                    autoColorHeight = FALSE, colorHeight = 0.5,cex.axis=.6,
                    cex.colorLabels = .7, abHeight = 0.95);
par(mar = c(8, 7.5, 1.5, .7));
labeledHeatmap(Matrix = pdat$numMat,
               xLabels = pdat$xLabels,xSymbols = pdat$xSymbols,
               yLabels = pdat$yLabels,ySymbols = pdat$ySymbols,
               colorLabels = TRUE, colors = blueWhiteRed(100)[50:100],
               textMatrix = pdat$textMat, cex.text = cex.textmatrix,
               setStdMargins = FALSE, cex.lab = cex.matrixlab,
               xColorWidth = 0.01,xLabelsAngle = 90,yColorWidth = 0.01,
               verticalSeparator.col = "lightgrey",
               verticalSeparator.interval = 1,
               horizontalSeparator.interval = 1,
               horizontalSeparator.col = "lightgrey",
               main = paste0("C. ",pdat$setname1," modules (rows) vs. ", pdat$setname2, " modules (columns)"),
               cex.main = cex.main);dev.off()

save(colorList,dendroList,
     modtab1,modtab2,
     meList,multiExpr,overlapList,
     scaledIMconlist,setPairs,goodIMC,statsObs,summaryTable,cex.main,
     cex.matrixlab,cex.textmatrix,
     keep,setLabels,file = "mp_overlap_strSpvC.rda")

