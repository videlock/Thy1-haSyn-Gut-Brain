
# example for colon (one month) shown. Uses processed data from step 1
# resulting top tables are in FinalProcessedData

library(Biobase)
library(limma)


setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))

grp="dc1m"

eset<-readRDS(file.path("../C_QCandNormalization",grp,"data/CorrectedEset.rds"))
GT<-factor(eset$GT)
design<-model.matrix(~0+GT)
colnames(design)<-levels(GT)
# head(design)
contrasts<-makeContrasts(Hem - WT, levels=design)
# contrasts

fit <- lmFit(eset, design)
fit.cont<-contrasts.fit(fit, contrasts)
fit.cont<-eBayes(fit.cont)
tt<-topTable(fit.cont,coef = "Hem - WT",adjust.method = "BH",number = "inf")
