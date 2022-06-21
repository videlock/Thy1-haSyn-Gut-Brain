rm(list = ls())
options(stringsAsFactors = FALSE);
library(WGCNA);
allowWGCNAThreads() 
library(flashClust);
library(gplots)
library(RColorBrewer)

setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))

load("modules.rda")

# Cluster module eigengenes---------
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");

pdf(file = "MEclustering.pdf",height=3.5,width=10
);par(mar=c(1,4,1,1)
);plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "",cex=.75);dev.off()



pdf("EigenNetwork.pdf")
par(mar=c(7,7,1,1))
labeledHeatmap(
  Matrix = cor(MEs),
  xLabels = names(MEs),
  yLabels = names(MEs),
  xSymbols =names(MEs),
  ySymbols = names(MEs),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  setStdMargins = FALSE,
  # xLabelsAngle = 0,
  xLabelsAdj = .8,
  cex.lab = .7
)
dev.off()



# ME traits cor plot------
source("../mylabeledHeatmap.R")

aso<-metaData$ASO==1
metaData$Age_ASO<-metaData$Age
metaData$Age_ASO[!aso]<-NA
metaData$Age_WT<-metaData$Age
metaData$Age_WT[aso]<-NA
rm(aso)

datTraitCols<-c("ASO_1m","ASO_3m","Age_ASO","Age_WT")
datTraits <- metaData[rownames(datExpr),datTraitCols]

names(datTraits)<-c("Thy1-haSyn v WT (1m)","Thy1-haSyn v WT (3m)", "3m v 1m (Thy1-haSyn)","3m v 1m (WT)")

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

moduleTraitPvalueBH<-apply(moduleTraitPvalue,2,function(x) p.adjust(x,method = "BH"))


textMatrixBH =  paste(signif(moduleTraitCor, 2),
                      " (",
                      signif(moduleTraitPvalueBH, 1),
                      ")",
                      sep = "")

dim(textMatrixBH) = dim(moduleTraitCor)

h.MEhm=6
w.MEhm=5
mar.MEhm= c(4, 6, 1, 0.5)
oma.MEhm=c(0,0,0,0)

pdf("MEtraitsCor.pdf",height=h.MEhm,width=w.MEhm
);par(mar = mar.MEhm,oma=oma.MEhm
);mylabeledHeatmap(
  Matrix = moduleTraitCor,
  textMatrix = textMatrixBH,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = gsub("ME","",names(MEs)),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  setStdMargins = FALSE,
  cex.text = 0.5,
  cex.lab.y = .6,
  cex.lab.x = 0.5,
  cex.legend = 0.5,
  xLabelsAngle = 45,
  xLabelsAdj = 0.9
);dev.off()



# modtab------
colnames(moduleTraitCor)<-paste(colnames(moduleTraitCor),"r",sep = "-")
colnames(moduleTraitPvalue)<-paste(colnames(moduleTraitPvalue),"p",sep = "-")
mods<-substring(colnames(MEs),3)
modsizes<-as.data.frame(table(factor(moduleColors, levels = mods)))
names(modsizes)<-c("moduleColor","Size")
cortab<-NULL

for(i in 1:ncol(moduleTraitCor)){
  cortab<-cbind(cortab,moduleTraitCor[,i],moduleTraitPvalueBH[,i])
}
cortab<-as.data.frame(cortab)
names(cortab)[seq(1,ncol(cortab)-1,2)]<-colnames(moduleTraitCor)
names(cortab)[seq(2,ncol(cortab),2)]<-colnames(moduleTraitPvalue)
modtab<-cbind(modsizes,cortab)

saveRDS(modtab,file = "modtab.rds")
# write.csv(modtab,"modtab.csv",row.names = F)

rm(MEDiss,METree,moduleTraitCor,moduleTraitPvalue,
   moduleTraitPvalueBH,textMatrixBH,modsizes,cortab)

# geneinfo------

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


geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
names(geneModuleMembership) = substring(names(MEs),3)
geneInfo<-makeGeneInfo(an[colnames(datExpr),],moduleColors,geneModuleMembership)


saveRDS(geneInfo,file = "geneINFOall.rds")
geneInfoNoGrey<-geneInfo[!geneInfo$moduleColor%in%"grey",] # no grey module

saveRDS(geneInfo[geneInfo$MM>=0.6,],file = "geneINFO_MM6.rds")

# write.csv(geneInfo[geneInfo$MM>=0.6,],"geneINFO.csv")

# GO ---

# create the run GO script
ntasks=3*4*nrow(modtab)

sink("runGO_str.R")
cat("#!/usr/bin/env Rscript",
    "#$ -cwd",
    paste0("#$ -o logs/go_$JOB_ID_$TASK_ID.txt"),
    "#$ -j y",
    paste0("#$ -t 1:",ntasks),sep = "\n")
sink()

system2(command = "cat",
        args = c(
          "runGO.R",
          ">>runGO_str.R"
        ))

# merge go results after running on cluster ----

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

