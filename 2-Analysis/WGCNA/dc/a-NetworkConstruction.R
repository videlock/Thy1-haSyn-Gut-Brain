rm(list = ls())

library(Biobase)
options(stringsAsFactors = FALSE);
library(WGCNA);
allowWGCNAThreads() 
library(flashClust);
library(gplots)
library(RColorBrewer)

setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))
eset<-readRDS("../../../1-RNAseqWorkflow/C_QCandNormalization/dc/data/CorrectedEset.rds") # path to corrected eset
an<-fData(eset)
an$Gene.description<-gsub("[[:space:]]\\[.+","",an$Gene.description)
names(an)[names(an)=="Chromosome.scaffold.name"]<-"Chr"
an<-an[!an$Chr%in%"X",]
datExpr<-exprs(eset)
keep=rownames(datExpr)%in%rownames(an)
datExpr<-datExpr[keep,]
datExpr<-as.data.frame(t(datExpr))
pheno<-pData(eset)
rownames(pheno)<-paste0("m",pheno$Mouse.ID)
rownames(datExpr)<-paste0("m",pheno$Mouse.ID)
pheno$ASO<-ifelse(pheno$GT=="Hem",1,0)
metaData<-pheno
t1m<-metaData$Time==1
metaData$ASO_1m<-as.numeric(as.factor(metaData$GT))-1
metaData$ASO_1m[!t1m]<-NA
metaData$ASO_3m<-as.numeric(as.factor(metaData$GT))-1
metaData$ASO_3m[t1m]<-NA
metaData$Age<-as.numeric(as.factor(metaData$Time))-1
metaData$GTAge<-metaData$GTtime

save(datExpr, metaData, an, file = "inputData_Raw.rda")

rm(list = ls())

# load data------
load("inputData_Raw.rda")

dataMatrix<-datExpr
dim(dataMatrix)

dataMatrix <- dataMatrix[,apply(dataMatrix,2,function(x) max(x[which(!is.na(x))])>0)]
dim(dataMatrix)



# if outliers were not checked before, run gsg and outliers here
gsg = goodSamplesGenes(dataMatrix, minNGenes=1, verbose = 3);#allok
gsg$allOK
dataMatrix_removeOutlier = dataMatrix[gsg$goodSamples, gsg$goodGenes]
dim(dataMatrix_removeOutlier)

metaData = metaData[rownames(metaData) %in% rownames(dataMatrix_removeOutlier),]
dim(metaData)

sdout <- 2
normadj <- (0.5+0.5*bicor(t(dataMatrix_removeOutlier), 
                          use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
sum(outliers) #0

dataMatrix_removeOutlier_keep <- dataMatrix_removeOutlier[!outliers,]
metaData = metaData[!outliers,]
dim(dataMatrix_removeOutlier_keep)


#clustering -------
datTraitCols<-c("ASO","ASO_1m","ASO_3m","Age")

sampleTree2 = hclust(dist(dataMatrix_removeOutlier_keep), method = "average")
datTraits <- metaData[rownames(dataMatrix_removeOutlier_keep),datTraitCols]
names(datTraits)<-c("Thy1-haSyn (all)","Thy1-haSyn (1m)","Thy1-haSyn (3m)","Age (3v1m)")

datTraits_colors<-datTraits
for(i in names(datTraits)){
  datTraits_colors[,i]<-numbers2colors(as.numeric(as.factor(datTraits[,i])))
}

colorPalette<-ifelse(ncol(datTraits_colors)<3,
                     brewer.pal(3, "Set1"),
                     brewer.pal(ncol(datTraits_colors), "Set1"))

names(colorPalette) = c(1:length(colorPalette))

pdf("sampleDendrogram.pdf",height=5,width=10)
plotDendroAndColors(sampleTree2, datTraits_colors,
                    groupLabels = names(datTraits_colors), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# save processed data-----
dir.create("processed_data")
datExpr <- dataMatrix_removeOutlier_keep
datTraits<-datTraits[rownames(metaData),]
# dim(datExpr)
save(datExpr, metaData, datTraits,an, file=file.path("processed_data","inputData.rda"))


# soft power thresholding-------
maxPout=5

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sftmax.file.name<-file.path("processed_data",paste0("signed.sft.thresh",maxPout,".rda"))


if(!file.exists(sftmax.file.name)){
  sft = pickSoftThreshold(datExpr, networkType="signed",corFnc="bicor",
                          corOptions = list(maxPOutliers =maxPout/100), 
                          powerVector = powers, verbose = 5)
  
  save(sft,file=sftmax.file.name)
}else{load(sftmax.file.name)}

pdf("softpower.pdf")

par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence, maxP outliers=",maxPout,"%"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")


# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity",maxPout,"%"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

rm(dataMatrix,dataMatrix_removeOutlier,dataMatrix_removeOutlier_keep,
   gsg,datTraits_colors,netsummary,normadj,outliers,
   colorPalette,powers,sftmax.file.name,sampleTree2,Z.K,Z.C,K,sdout,
   sft,C,cex1)

# dendo with cutting parameters---------
sp=8

adjacency = adjacency(datExpr, power = sp, type = "signed",
                      corFnc="bicor",corOptions = list(maxPOutliers =maxPout/100))
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = flashClust(as.dist(dissTOM), method = "average");


traitmat<-datTraits
geneSigs=matrix(NA,nrow=length(datTraitCols),ncol=ncol(datExpr)) # create a vector to hold the data
colnames(geneSigs)=colnames(datExpr);rownames(geneSigs)=names(datTraits)

for(i in 1:ncol(geneSigs)) {
  # calculate adjusted R^2s square-root for categorical variables (factors)
  exprvec=as.numeric(datExpr[,i]) # get the expression vector for ith gene
  for(j in names(datTraits)){
    geneSigs[j,i]<-sqrt(max(
      summary(lm(exprvec~as.factor(traitmat[,j])))$adj.r.squared,0))
  }
}

# convert to colors
geneSigsColor=matrix(NA,nrow=nrow(geneSigs),ncol=ncol(datExpr)) # create a vector to hold the data
rownames(geneSigsColor)<-rownames(geneSigs)
colnames(geneSigsColor)<-colnames(datExpr)

for(i in names(datTraits)){
  geneSigsColor[i,]<-numbers2colors(as.numeric(geneSigs[i,]), 
                                    signed=FALSE, centered=FALSE, 
                                    blueWhiteRed(100)[51:100], lim=c(0,1))
}


saveRDS(geneTree,file ="geneTree.rds")
saveRDS(geneSigsColor,file ="geneSigsColor.rds")

MMSs<-c(30,60)
DSs<-c(1,2)
DTs<-c(0.1,0.2,0.25)
cutHeights<-c(0.995)

mColorh <- mLabelh <- colorLabels <- NULL

for(ch0 in cutHeights){
  for (mms0 in MMSs) {
    for (ds0 in DSs) {
      for (dthresh0 in DTs) {
        ch.lab=gsub("0\\.","",ch0)
        dthresh.lab=gsub("0\\.","",dthresh0)
        paramName<-paste0("ch",ch.lab,"ds",ds0,"mm",mms0,"dt",dthresh.lab)
        dynamicMods = cutreeDynamic(cutHeight = ch0,dendro = geneTree,
                                    minClusterSize = mms0,
                                    deepSplit = ds0,
                                    distM = as.matrix(dissTOM))
        
        dynamicColors = labels2colors(dynamicMods)
        
        merged = mergeCloseModules(datExpr, dynamicColors,
                                   cutHeight = dthresh0,
                                   verbose = 0)
        
        moduleColors = merged$colors
        nmods=length(merged$dendro$labels)
        
        mColorh0 <- moduleColors
        
        mColorh <- cbind(mColorh,mColorh0)
        
        mLabelh0 <- paste0("ch=",ch.lab," ds=",ds0," mms=",mms0,
                           "\ndt=",dthresh.lab," nmods=",nmods)
        mLabelh <- c(mLabelh,mLabelh0)
        
      }
    }
  }
}


mColorh1=cbind(mColorh,t(geneSigsColor))
mLabelh1=c(mLabelh,rownames(geneSigsColor))

pdf("Signed_Dendro_parameters.pdf",height = 10,width = 20)
plotDendroAndColors(geneTree,mColorh1,groupLabels=mLabelh1,autoColorHeight = F,
                    addGuide=TRUE, dendroLabels=FALSE,colorHeight = .65,
                    main="Dendrogram With Different Module Cutting Parameters")
dev.off()

rm(adjacency,TOM,ch.lab,ch0,colorLabels,cutHeights,ds0,DSs,dthresh.lab,
   dthresh0,DTs,exprvec,i,j,dynamicMods,dynamicColors,mColorh0,mColorh,
   mColorh1,mLabelh0,mLabelh,mLabelh1,mms0,MMSs,moduleColors,nmods,paramName)


# use selected parameters--------
ch<-0.995
mms<-60
ds<-1
dthresh<-0.25

dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight = ch,
                            minClusterSize = mms, 
                            deepSplit = ds, distM = as.matrix(dissTOM))
dynamicColors = labels2colors(dynamicMods)
merged = mergeCloseModules(datExpr, dynamicColors, cutHeight = dthresh, verbose = 0)
moduleColors = merged$colors


mColorh <- cbind(moduleColors,t(geneSigsColor))
mLabelh <- c("Merged Colors",rownames(geneSigsColor))

pdf("Signed_Dendro_params_Final.pdf",height=10,width=16)
plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh, 
                    addGuide=TRUE,dendroLabels=FALSE, 
                    main= paste("Signed bicor network with power = ",sp,"mms=",mms,"ds=",ds,"dthresh=",dthresh));

dev.off()


MEs0 = moduleEigengenes(datExpr, moduleColors,excludeGrey = F)$eigengenes
MEs = orderMEs(MEs0)
rownames(MEs) = rownames(datExpr)
names(moduleColors) <- colnames(datExpr)
save(geneTree,moduleColors,MEs,sp,mms,ch,dthresh,maxPout,
     sd,datTraitCols,metaData,datExpr,an,file="modules.rda")

# to have modules match paper modules, run the following
paperColors<-c('darkred', 'violet', 'blue', 'brown', 'lightgreen', 'grey60', 'saddlebrown', 'red', 'lightcyan', 'purple', 'salmon', 'orange', 'sienna3', 'lightyellow', 'yellowgreen', 'tan', 'white', 'black', 'darkorange', 'greenyellow', 'royalblue', 'darkmagenta', 'skyblue', 'darkolivegreen', 'magenta', 'steelblue', 'midnightblue', 'paleturquoise', 'yellow', 'darkturquoise', 'turquoise', 'green', 'cyan', 'darkgreen', 'darkgrey', 'pink')
mods<-substring(colnames(MEs),3)

modConvert<-data.frame(old=mods,new=paperColors)
moduleColors.new<-character(length = length(moduleColors))

for(i in 1:nrow(modConvert)){
  moduleColors.new[moduleColors==modConvert$old[i]]<-modConvert$new[i]
}

MEs0 = moduleEigengenes(datExpr, moduleColors.new,excludeGrey = F)$eigengenes
MEs = orderMEs(MEs0)
rownames(MEs) = rownames(datExpr)
names(moduleColors) <- colnames(datExpr)
save(geneTree,moduleColors,MEs,sp,mms,ch,dthresh,maxPout,
     sd,datTraitCols,metaData,datExpr,an,file="modules.rda")
