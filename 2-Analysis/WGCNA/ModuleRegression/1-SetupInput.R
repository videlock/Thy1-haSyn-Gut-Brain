rm(list = ls())
options(stringsAsFactors = FALSE);
library(WGCNA);
allowWGCNAThreads()
library(flashClust);
library(gplots)
library(RColorBrewer)
library(limma)
library(ggplot2)
library(MASS)

setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))


# setup-----------
load("../Colon/modules.rda")
load("../Colon/processed_data/inputData.rda")
names(MEs)<-gsub("ME","dc.",names(MEs))
dc.MEs<-MEs
dc.dat<-cbind(metaData[,c(3:4,14:18)],dc.MEs)
dc.dat$Genotype<-ifelse(dc.dat$GT=="Hem","Thy1-haSyn","WT")
dc.dat$Timepoint<-ifelse(dc.dat$Time==1,"1 month","3 months")

load("../Striatum/modules.rda")
load("../Striatum/processed_data/inputData.rda")

names(MEs)<-gsub("ME","str.",names(MEs))
str.MEs<-MEs
str.dat<-cbind(metaData[,c(3:4,13:17)],MEs)
str.dat$Genotype<-ifelse(str.dat$GT=="Hem","Thy1-haSyn","WT")
str.dat$Timepoint<-ifelse(str.dat$Time==1,"1 month","3 months")


commonsamps<-intersect(rownames(dc.dat),rownames(str.dat))

bg.dat<-cbind(dc.dat[commonsamps,],str.MEs[commonsamps,])

table(bg.dat$GTAge)
# Hem_1 Hem_3  WT_1  WT_3 
# 4     6     3     3 

save(dc.dat,str.dat,bg.dat,str.MEs,dc.MEs, file = "regInputData.rda")

