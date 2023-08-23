
rm(list=ls());

options(echo=TRUE)
options(stringsAsFactors = FALSE);
library(WGCNA);
allowWGCNAThreads() 
library(flashClust);
library(gplots)
library(RColorBrewer)

setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))

# load networks or start with input below
load("../Colon/processed_data/inputData.rda")

datexpr.dc<-datExpr

datTraitNumericCols<-c("ASO","Age","ASO_1m","ASO_3m")
datTraits.dc <- metaData[rownames(datExpr),datTraitNumericCols]

an.dc<-an
metaData.dc<-metaData

load("../Striatum/processed_data/inputData.rda")

datexpr.str<-datExpr
datTraits.str <- metaData[rownames(datExpr),datTraitNumericCols]
an.str<-an
metaData.str<-metaData

datTraitsList<-list(dc=datTraits.dc,str=datTraits.str)
metadataList<-list(dc=metaData.dc,str=metaData.str)

Genes<-intersect(names(datexpr.dc),names(datexpr.str))
datexpr.dc<-datexpr.dc[,Genes]
datexpr.str<-datexpr.str[,Genes]

multiExpr<-multiData(dc=datexpr.dc, str=datexpr.str)

an<-an[Genes,]

save(multiExpr,an,datTraitsList,metadataList,file = "inputData.rda")

rm(list = ls())

load("inputData.rda")

net = blockwiseConsensusModules(
  multiExpr, 
  power = c(8,12), 
  corType = "bicor",
  saveIndividualTOMs = F,
  saveConsensusTOMs = F,
  minModuleSize = 60, 
  deepSplit = 1,
  detectCutHeight = .995,
  maxPOutliers = c(0.05,0.01),
  pamRespectsDendro = FALSE, 
  networkType = "signed",
  mergeCutHeight = 0.25, 
  numericLabels = FALSE,
  minKMEtoStay = 0.2, 
  verbose = 5,
  maxBlockSize = 50000)


setLabels = c("Colon", "Striatum")
shortLabels = c("dc", "str")
mods<-substring(names(net$multiMEs$dc$data),3)
consMEs.ConsOrder<-consensusOrderMEs(net$multiMEs)

save(multiExpr,net,
     setLabels,shortLabels,consMEs.ConsOrder,mods,
     file="modules.rda")

mar.eigennet=c(1.5,7.5,0,0.5)
oma.eigennet=c(0,0,1,0)
mar.MEhm= c(1.25, 7.5, 0.25, 0.5)
oma.MEhm=c(0,0,1,0)
h.MEhm=7
w.MEhm=6
mar.dendro=c(0,6,4,0)
mar.net=c(2,5,1,0)

source("../mylabeledHeatmap.R")
source("../myplotEigengenNetworks.R")

#eigenplots----
pdf("consEigenPlots.pdf"
);myplotEigengeneNetworks(
  consMEs.ConsOrder, 
  setLabels, 
  letterSubPlots = FALSE, Letters = NULL, 
  excludeGrey = TRUE, greyLabel = "grey", 
  plotDendrograms = FALSE, plotHeatmaps = TRUE, 
  setMargins = TRUE, marDendro = NULL,
  colorLabels = TRUE, signed = TRUE, 
  heatmapColors = NULL, 
  plotAdjacency = FALSE,
  printAdjacency = FALSE, cex.adjacency = 0.9,
  coloredBarplot = TRUE, barplotMeans = FALSE, barplotErrors = FALSE, 
  plotPreservation = "standard",
  printPreservation = FALSE, cex.preservation = 0.9,
  marHeatmap = c(3,3,2,1),
  zlimPreservation = c(0, 1),
  xLabelsAngle = 90
);dev.off()


# MEcorplots---------
datTraitNumericCols<-c("ASO","ASO_1m","ASO_3m","Age")
dc.datTraits <- datTraitsList$dc
dc.datTraits <- dc.datTraits[,datTraitNumericCols]
names(dc.datTraits)<-c("Thy1-haSyn (all)","Thy1-haSyn (1m)","Thy1-haSyn (3m)","Age (3v1m)")
str.datTraits <- datTraitsList$str
str.datTraits <- str.datTraits[,datTraitNumericCols]
names(str.datTraits)<-c("Thy1-haSyn (all)","Thy1-haSyn (1m)","Thy1-haSyn (3m)","Age (3v1m)")

dc.datExpr=multiExpr[["dc"]][["data"]]
dc.nGenes = ncol(dc.datExpr);dc.nSamples = nrow(dc.datExpr)
str.datExpr=multiExpr[["str"]][["data"]]
str.nGenes = ncol(str.datExpr);str.nSamples = nrow(str.datExpr)

dc.MEs=consMEs.ConsOrder[["dc"]][["data"]]
dc.MEs<-dc.MEs[,!colnames(dc.MEs)%in%"MEgrey"]

str.MEs=consMEs.ConsOrder[["str"]][["data"]]
str.MEs<-str.MEs[,!colnames(str.MEs)%in%"MEgrey"]

dc.moduleTraitCor = cor(dc.MEs, dc.datTraits, use = "p")
dc.moduleTraitPvalue = corPvalueStudent(dc.moduleTraitCor, dc.nSamples)
dc.textMatrix =  paste(signif(dc.moduleTraitCor, 2),
                       " (",
                       signif(dc.moduleTraitPvalue, 1),
                       ")",
                       sep = "")
dim(dc.textMatrix) = dim(dc.moduleTraitCor)
dc.main<-"Colon"

str.moduleTraitCor = cor(str.MEs, str.datTraits, use = "p")
str.moduleTraitPvalue = corPvalueStudent(str.moduleTraitCor, str.nSamples)
str.textMatrix =  paste(signif(str.moduleTraitCor, 2),
                        " (",
                        signif(str.moduleTraitPvalue, 1),
                        ")",
                        sep = "")
dim(str.textMatrix) = dim(str.moduleTraitCor)
str.main<-"Striatum"


pdf("consMEtraitsCor.pdf",height=4.5,width=10
);layout(mat = matrix(c(1,2),nrow = 1),widths = c(2, 2)
);par(mar = c(1.5, 7, 1.5, 0)
);mylabeledHeatmap(
  Matrix = dc.moduleTraitCor, textMatrix = dc.textMatrix,
  xLabels = names(dc.datTraits), yLabels = names(dc.MEs),
  ySymbols = substring(names(dc.MEs),3),
  colorLabels = TRUE,plotLegend = F,
  colors = blueWhiteRed(50),
  setStdMargins = FALSE, cex.text = .75,cex.legend = 0.7,
  cex.lab.y = 1, cex.lab.x = 0.7, xLabelsAngle = 0, xLabelsAdj = 0.5,
  cex.main=0.75,keepLegendSpace = F,
);title(dc.main,adj=0.44
);par(mar = c(1.5, 2, 1.5, 0.5)
);mylabeledHeatmap(
  Matrix = str.moduleTraitCor, textMatrix = str.textMatrix,
  xLabels = names(str.datTraits), yLabels = names(str.MEs),
  # ySymbols = substring(names(str.MEs),3),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  setStdMargins = FALSE, cex.text = .75,cex.legend = 0.7,
  cex.lab.y = 1, cex.lab.x = 0.7, xLabelsAngle = 0, xLabelsAdj = 0.5,
  cex.main=0.75
);title(str.main,adj=0.45
);dev.off()

# combined eigen and ME----------

layoutmat=rbind(c(1,3),
                c(2,4),
                c(5,7),
                c(6,8),
                c(9,0),
                c(10,11))

cex.subtitle=1.5
adj.subtitle=c(0,0.5)
x.subtitle=0;y.subtitle=0.5
cex.lab.x = 0.75

pdf("consensusNet.pdf",height=8,width=7.5);

layout(mat = layoutmat,
       widths = c(2, 2),heights = c(0.4,4,0.4,4,0.4,5));
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(x.subtitle,y.subtitle,"A. Colon",cex = cex.subtitle, font = 2,adj = adj.subtitle)

tadj=-.1
cex.adjacency = 0.9
plotPreservation = "standard"
zlimPreservation = c(0,1)
cex.preservation = 0.9
greyLabel = "grey"
multiME<-consMEs.ConsOrder
heatmapColors = blueWhiteRed(50)
nSets = length(multiME)


for (set in 1:nSets){
  multiME[[set]]$data = 
    multiME[[set]]$data[ , substring(names(multiME[[set]]$data),3)!=greyLabel]
}
labels = names(multiME[[set]]$data);
uselabels = labels[substring(labels,3)!=greyLabel];
corME = cor(multiME[[set]]$data[substring(labels,3)!=greyLabel,
                                substring(labels,3)!=greyLabel], use="p");
disME = as.dist(1-corME);
clust = fastcluster::hclust(disME, method = "average");  
plotLabels = uselabels;

# colon
par(mar = c(2,2,0,0))
i.row=1;i.col=1
nModules = dim(multiME[[i.col]]$data)[2]
textMat = NULL;
corME = cor(multiME[[i.col]]$data, use="p") 
pME = corPvalueFisher(corME, nrow(multiME[[i.col]]$data))
mylabeledHeatmap((1+corME)/2, names(multiME[[i.col]]$data), names(multiME[[i.col]]$data),
                 # main=paste(letter, setLabels[[i.col]]), 
                 invertColors=FALSE,
                 zlim=c(0,1.0),
                 colorLabels = T, colors = heatmapColors, 
                 setStdMargins = FALSE, cex.legend = 0.7,
                 textMatrix = textMat, cex.text = cex.adjacency);

# barplot
i.row=1;i.col=2
corME1 = cor(multiME[[i.col]]$data, use="p");
corME2 = cor(multiME[[i.row]]$data, use="p");
cor.dif = (corME1 - corME2)/2;
d = tanh((corME1 - corME2) / (abs(corME1) + abs(corME2))^2);
dispd = cor.dif
dp = 1-abs(cor.dif); 
method = "Preservation:";
diag(dp) = 0
sum_dp = sum(dp[upper.tri(dp)]);

par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(x.subtitle,y.subtitle,labels = paste("B.", method, "sum = ", signif(sum_dp,3)),cex = cex.subtitle, font = 2,adj = adj.subtitle)

par(mar = c(2,2,0,0))
labeledBarplot(dp, names(multiME[[i.col]]$data),
               # main=paste(letter, method, "sum = ", signif(sum_dp,3)), 
               ylim=c(0,dim(dp)[[1]]),cex.lab = 0.5,
               colorLabels = T, colored = T, 
               setStdMargins = FALSE)

# preservation plot
i.row=2;i.col=1
corME1 = cor(multiME[[i.col]]$data, use="p");
corME2 = cor(multiME[[i.row]]$data, use="p");
cor.dif = (corME1 - corME2)/2;
d = tanh((corME1 - corME2) / (abs(corME1) + abs(corME2))^2);
dispd = cor.dif
half = as.integer(length(heatmapColors)/2);
range = c(half:length(heatmapColors)); 
halfColors = heatmapColors[range];
printMtx = NULL

par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(x.subtitle,y.subtitle,"C. Preservation",cex = cex.subtitle, font = 2,adj = adj.subtitle)

par(mar = c(2,2,0,0))
mylabeledHeatmap(1-abs(dispd), names(multiME[[i.col]]$data), names(multiME[[i.col]]$data), 
                 # main = main, 
                 invertColors=FALSE,
                 colorLabels = T, zlim = zlimPreservation, colors = halfColors,
                 setStdMargins = FALSE, cex.legend = 0.7,
                 textMatrix = printMtx, cex.text = cex.preservation);

# striatum
i.row=2;i.col=2

par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(x.subtitle,y.subtitle,"D. Striatum",cex = cex.subtitle, font = 2,adj = adj.subtitle)

par(mar = c(2,2,0,0))

nModules = dim(multiME[[i.col]]$data)[2]
textMat = NULL;
corME = cor(multiME[[i.col]]$data, use="p") 
pME = corPvalueFisher(corME, nrow(multiME[[i.col]]$data))
mylabeledHeatmap((1+corME)/2, names(multiME[[i.col]]$data), names(multiME[[i.col]]$data),
                 # main=paste(letter, setLabels[[i.col]]), 
                 invertColors=FALSE,
                 zlim=c(0,1.0),
                 colorLabels = T, colors = heatmapColors, 
                 setStdMargins = FALSE, cex.legend = 0.7,
                 textMatrix = textMat, cex.text = cex.adjacency);


par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(x.subtitle,y.subtitle,"E. ME trait correlations",cex = cex.subtitle, font = 2,adj = adj.subtitle)

par(mar = c(1.5, 7, 2, 0));
mylabeledHeatmap(
  Matrix = dc.moduleTraitCor, textMatrix = dc.textMatrix,
  xLabels = names(dc.datTraits), yLabels = names(dc.MEs),
  ySymbols = substring(names(dc.MEs),3),
  colorLabels = TRUE,plotLegend = F,
  colors = blueWhiteRed(50),
  setStdMargins = FALSE, cex.text = .75,cex.legend = 0.7,
  cex.lab.y = 1, xLabelsAngle = 0, xLabelsAdj = 0.5,
  cex.main=0.75,keepLegendSpace = F, cex.lab.x = cex.lab.x
);title(dc.main,adj=0.44
);

par(mar = c(1.5, 2, 2, 0.5)
);mylabeledHeatmap(
  Matrix = str.moduleTraitCor, textMatrix = str.textMatrix,
  xLabels = names(str.datTraits), yLabels = names(str.MEs),
  # ySymbols = substring(names(str.MEs),3),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  setStdMargins = FALSE, cex.text = .75,cex.legend = 0.7,
  cex.lab.y = 1, xLabelsAngle = 0, xLabelsAdj = 0.5,
  cex.main=0.75, cex.lab.x = cex.lab.x
);title(str.main,adj=0.45
);dev.off()



makeGeneInfo<-function(geneData,modules,geneModMemb){
  geneInfo0 = data.frame(ENSG = rownames(geneData),
                         geneSymbol = geneData$Gene.name,
                         EntrezID = geneData$NCBI.gene.ID,
                         Name=geneData$Gene.description,
                         moduleColor = modules)
  Modmemb<-numeric()
  for(i in 1:nrow(geneInfo0)){
    modColor<-modules[i]
    Modmemb[i]<-geneModMemb[i,modColor]
  }
  geneInfo0$MM<-Modmemb
  geneOrder = order(factor(geneInfo0$moduleColor,levels = names(geneModMemb)), -geneInfo0$MM)
  
  geneInfo<-geneInfo0[geneOrder,]
  geneInfo<-geneInfo[geneInfo$geneSymbol!="",]
  return(geneInfo)
}


geneModuleMembership = as.data.frame(cor(multiExpr$dc$data,
                                         consMEs.ConsOrder$dc$data, use = "p"))
names(geneModuleMembership) = substring(names(consMEs.ConsOrder$dc$data),3)
geneInfo.dc<-makeGeneInfo(an[colnames(multiExpr$dc$data),],
                          net$colors,geneModuleMembership)


geneModuleMembership = as.data.frame(cor(multiExpr$str$data,
                                         consMEs.ConsOrder$str$data, use = "p"))
names(geneModuleMembership) = substring(names(consMEs.ConsOrder$str$data),3)
geneInfo.str<-makeGeneInfo(an[colnames(multiExpr$str$data),],
                           net$colors,geneModuleMembership)

# both MM columns
geneInfo.both<-geneInfo.dc[order(geneInfo.dc$ENSG),]

names(geneInfo.both)[6]<-"MM.dc"
geneInfo.str2<-geneInfo.str[order(geneInfo.str$ENSG),]
geneInfo.both$MM.str<-geneInfo.str2$MM
geneInfo.both$avgMM<-rowMeans(geneInfo.both[,c("MM.dc","MM.str")])


geneOrder = order(factor(geneInfo.both$moduleColor,levels = names(geneModuleMembership)), -geneInfo.both$avgMM)
geneInfo.both<-geneInfo.both[geneOrder,]

geneInfo.both_all<-geneInfo.both
geneInfo.both<-geneInfo.both[(geneInfo.both$MM.dc>=0.6)|(geneInfo.both$MM.str>=0.6),]

write.csv(geneInfo.both,"geneInfoboth.csv",row.names = F)
write.csv(geneInfo.both_all,"geneInfobothAll.csv",row.names = F)
saveRDS(geneInfo.both,file = "geneINFOboth_MM6.rds")
saveRDS(geneInfo.both_all,file = "geneINFOallgenes.rds")

# mod table-----

modtablist<-list()
for(tis in names(datTraitsList)){
  datTraits <- datTraitsList[[tis]]
  datExpr=multiExpr[[tis]][["data"]]
  nSamples = nrow(datExpr)
  MEs=consMEs.ConsOrder[[tis]][["data"]]
  MEs<-MEs[,!colnames(MEs)%in%"MEgrey"]
  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  colnames(moduleTraitCor)<-paste(names(datTraits),tis,sep = "-")
  colnames(moduleTraitCor)<-paste(colnames(moduleTraitCor),"r",sep = "-")
  colnames(moduleTraitPvalue)<-paste(names(datTraits),tis,sep = "-")
  colnames(moduleTraitPvalue)<-paste(colnames(moduleTraitPvalue),"p",sep = "-")
  mods<-substring(colnames(MEs),3)
  modsizes<-as.data.frame(table(factor(net$colors, levels = mods)))
  names(modsizes)<-c("moduleColor","Size")

  
  cortab<-NULL
  
  for(i in 1:ncol(moduleTraitCor)){
    cortab<-cbind(cortab,moduleTraitCor[,i],moduleTraitPvalue[,i])
  }
  cortab<-as.data.frame(cortab)
  names(cortab)[seq(1,ncol(cortab)-1,2)]<-colnames(moduleTraitCor)
  names(cortab)[seq(2,ncol(cortab),2)]<-colnames(moduleTraitPvalue)
  
  modtablist[[tis]]<-cbind(modsizes,cortab)
}
modtab<-cbind(modtablist[[1]][,c(1,2,3,4)],
                modtablist[[2]][,c(3,4)],
                modtablist[[1]][,c(5,6)],
                modtablist[[2]][,c(5,6)],
                modtablist[[1]][,c(7,8)],
                modtablist[[2]][,c(7,8)],
                modtablist[[1]][,c(9,10)],
                modtablist[[2]][,c(9,10)])


saveRDS(modtab,file = "modtab.rds")
write.csv(modtab,"modtab.csv",row.names = F)

# GO ---

ntasks=3*4*nrow(modtab)

sink("runGO_BGcons.R")
cat("#!/usr/bin/env Rscript",
    "#$ -cwd",
    paste0("#$ -o logs/go_$JOB_ID_$TASK_ID.txt"),
    "#$ -j y",
    paste0("#$ -t 1:",ntasks),sep = "\n")
sink()

system2(command = "cat",
        args = c(
          "rungo_cons.R",
          ">>runGO_BGcons.R"
        ))



# merge go results ----

onts<-c("BP","CC","MF")
MMthreshs<-c("6","7","8","all")
mods<-as.character(modtab$moduleColor)


param.df<-data.frame()

for(mod in mods){
  for(ont in onts){
    for(MMthresh in MMthreshs){
      param.df<-rbind(param.df,
                      c(mod,ont,MMthresh))
    }
  }
}


names(param.df)<-c("mod","ont","MMthresh")

param.df$modpath<-file.path("GO",
                            param.df$ont,
                            paste0("MM",param.df$MMthresh),
                            paste0(param.df$mod,".rds")
)
keep<-unlist(lapply(param.df$modpath,file.exists))
param.df<-param.df[keep,]


makeGOtab<-function(x){
  gotab<-NULL
  # golist<-lapply(x$modpath,readRDS)
  for(i in 1:nrow(x)){
    mod=x$mod[i]
    ont=x$ont[i]
    mm=paste0("MM",x$MMthresh[i])
    go<-readRDS(x$modpath[i])
    names(go)[1]<-"GOID"
    go$moduleColor<-mod
    go$ont<-ont
    go$mm<-mm
    go<-go[go$Size<1000,]
    go<-go[go$Count>3,]
    if(nrow(go)>0){
      gotab<-rbind(gotab,go)
    }
  }
  
  gotabMods<-split(gotab,gotab$moduleColor)
  
  gotab2<-NULL
  for(mod in names(gotabMods)){
    modtab<-gotabMods[[mod]]
    modtab$Term2<-paste(modtab$ont,modtab$Term,sep = ": ")
    modtab$termPval<-ifelse(modtab$Pvalue<0.005,
                            paste0(modtab$Term2, " (<0.005)"),
                            paste0(modtab$Term2, " (", round(modtab$Pvalue,digits = 3),")"))
    modtab<-modtab[order(modtab$Pvalue),]
    y<-modtab[modtab$Pvalue<0.005,]
    if(nrow(y)>5){
      moddat<-y
    }else if(nrow(x)>20){
      moddat<-x[1:20,]
    }else{moddat<-modtab}
    gotab2<-rbind(gotab2,modtab)
  }
  
  return(gotab2)
}


go<-makeGOtab(param.df)
go$moduleColor<-factor(go$moduleColor,levels=modtab$moduleColor)
go<-go[order(go$moduleColor),]
go$moduleColor<-as.character(go$moduleColor)
saveRDS(go,file="GOmulti.rds")
write.csv(go,row.names = F,file = "GOmulti.csv")


