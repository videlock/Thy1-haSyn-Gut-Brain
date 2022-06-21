rm(list = ls())

# the libraries needed and the actual code for this section is within convenience functions in QCfunctions.R 

source("QCfunctions.R") 
grp="dc1m"
outname<-grp


# set wd to QC if needed
# setwd(list.dirs(recursive = F)[grep("C_QC",list.dirs(recursive = F))])

# setup directory structure
dirSetup(out.name = outname)

# load eset and restrict to group
eset<-loadEset(Grp = grp, path.eset = "data/rawCountEset.rds")
dim(eset);table(eset$GT)


# filter, outliers, variance stabilization ---------------

# filter
filteredEset<-filterEset(Eset = eset,
                         CountsGrThan = 5,
                         inXpercent = 30,
                         savePath = file.path(datdir,"filteredEset.rds"))
dim(filteredEset)


# outliers

dat.exprs<-getDatExprs(filteredEset)
dataMatrix_removeOutliers<-getGSG(dat.exprs)


ZKoutliers<-findZKoutliers(DatExprs = dataMatrix_removeOutliers, sdout = 2)
plotSampleTree(DatExprs = dataMatrix_removeOutliers,
               savePath = file.path(figdir,"SampleTree.pdf"))

# repeat filtering removing the two outliers
eset<-loadEset(Grp = grp, 
               path.eset = "data/rawCountEset.rds", 
               Outliers = c("7648_dc","7675_dc"))

filteredEset<-filterEset(Eset = eset,
                         CountsGrThan = 5,
                         inXpercent = 30,
                         savePath = file.path(datdir,"filteredEset.rds"))
dim(filteredEset)


# variance stabilization ---------------

cds <- DESeqDataSetFromMatrix(countData = exprs(filteredEset),
                              colData = pData(filteredEset),
                              design = ~ GT)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
vsd <- getVarianceStabilizedData(cds)
normDat<-vsd
saveRDS(normDat,file=file.path(datdir,"NormCounts.rds"))

#make into eset
filtNormEset<-filteredEset
exprs(filtNormEset)<-normDat
filtNormEsetFile<-file.path(datdir,"filtNormEset.rds")
saveRDS(filtNormEset,file=filtNormEsetFile)

# density plots -------------------

printDensityAndMeanVarPlots(fig.dir = figdir,
                            rawData = exprs(eset),
                            FiltData = exprs(filteredEset),
                            FiltNormData = exprs(filtNormEset),
                            DispData = cds
                            )

# Removal of sequencing artifacts----------
conds0<-pData(filtNormEset)
conds0$sample<-rownames(conds0)
conds<-conds0[,c("sample","GT")]
target<-makeTargetDat(picard.path="data/picardTable.rds",Conds = conds)
normExpr<-as.data.frame(exprs(filtNormEset))

# sequencing PCs
target<-addSeqPCs(Target = target, Conds = conds,
                  norm.Expr = normExpr)

# expression PCs
PCAdatlist<-list()
PCAdatlist$seq$topvar<-seqtopvar
PCAdatlist$pre<-PCvar(normExpr)

plotvar(pclist = PCAdatlist, expPCs = "pre",
        outdir = figdir,fname = "screes.pdf")

# pairsdat for correlation matrix
cond=labels2colors(as.numeric(factor(target$GT))) 
n=dim(conds)[2]
pairsdat <- data.frame(
  GT = as.factor(target$GT),
  Seq_PC1 = as.numeric(target$Seq.PC1),
  Seq_PC2 = as.numeric(target$Seq.PC2),
  Seq_PC3 = as.numeric(target$Seq.PC3),
  Seq_PC4 = as.numeric(target$Seq.PC4),
  Seq_PC5 = as.numeric(target$Seq.PC5),
  Seq_PC6 = as.numeric(target$Seq.PC6),
  Seq_PC7 = as.numeric(target$Seq.PC7)
)


# pre-correction Plots ----------

mdslist<-list()
mdslist$pre<-makemds.set(data = normDat)

mdsplot1v2(mdslist$pre,figdir,"MDS_pre.pdf")
mdsplotAll(mdslist$pre,figdir,"MDS_pre_all.pdf")

plothms(pclist = PCAdatlist, expPCs = "pre",
        Conds = conds, Target = target,
        outdir = figdir, fname = "hms1.pdf")


cormatrix(pclist = PCAdatlist, expPCs = "pre", fname = "CorrPlot_Pre.pdf",
          Narrow = F)


# correction setup -----------------
GT<-as.numeric(as.factor(target$GT))-1
seqPC1 <- as.numeric(target$Seq.PC1)
seqPC2 <- as.numeric(target$Seq.PC2)
seqPC3 <- as.numeric(target$Seq.PC3)
seqPC4 <- as.numeric(target$Seq.PC4)
seqPC5 <- as.numeric(target$Seq.PC5)
seqPC6 <- as.numeric(target$Seq.PC6)
seqPC7 <- as.numeric(target$Seq.PC7)

dat.cor<-list()


# remove sequencing pcs 1, 2, & 3 -------------
modVars<-c("GT","seqPC1","seqPC2","seqPC3")
corname<-"SeqPC123"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsCompPlot(mds.list = mdslist[c("pre",corname)],fname = "MDS_Comp-Pre-SeqPC123.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC123_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorrSeqPC123.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_CorrSeqPC123.pdf")


# remove sequencing pcs 1, 2, 3 & 7 -------------
modVars<-c("GT","seqPC1","seqPC2","seqPC3","seqPC7")
corname<-"SeqPC1237"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsCompPlot(mds.list = mdslist[c("pre",corname)],fname = "MDS_Comp-Pre-SeqPC1237.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC1237_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorrSeqPC1237.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_CorrSeqPC1237.pdf")


# save corrected as eset------------
eset.cor<-filtNormEset
dim(pData(eset.cor))
dim(eset.cor)
dim(dat.cor$SeqPC1237)
exprs(eset.cor)<-as.matrix(dat.cor$SeqPC1237)
saveRDS(eset.cor, file = file.path(datdir,"CorrectedEset.rds"))
saveRDS(target, file = file.path(datdir,"target.rds"))
