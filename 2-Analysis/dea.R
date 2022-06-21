# setup ---------
library(Biobase)
library(limma)

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
