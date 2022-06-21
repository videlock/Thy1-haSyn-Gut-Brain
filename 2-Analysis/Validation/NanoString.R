
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))


## load data
library(NanoStringNorm)
load("NSraw.rda")
Nano_SampleInfo<-readRDS("NanoSampleInfo.rds")
ns.data<-NanoString.mRNA$x[,-c(1:3)]
ns.anno <- NanoString.mRNA$x[,c(1:3)]
ns.header<-NanoString.mRNA$header
traits<-data.frame(gt=ifelse(Nano_SampleInfo$GT=="Hem",1,2),
                   time=ifelse(Nano_SampleInfo$Time=="4 wks",1,2))
ns.anno$Name[ns.anno$Name=="Cacna2s3"]<-"Cacna2d3"


## Normalize just code count and background (positive, negative controls)
nsNorm.cc.bg<-NanoStringNorm(ns.data,
                             anno = ns.anno,
                             header = ns.header,
                             CodeCount = "geo.mean",
                             Background = "mean.2sd",
                             round.values = T,take.log = T, traits = traits)



### genes with missing in over 90%
endog_to_rem<-rownames(
  nsNorm.cc.bg$gene.summary.stats.norm[
    nsNorm.cc.bg$gene.summary.stats.norm$Missing>90,
  ])
endog_to_rem<-endog_to_rem[endog_to_rem%in%endog]
endog_to_rem

### hks to keep - low CV, range of expression level
hksToKeep<-c("Actn1","Gusb","Hprt","Rplp0")

exclude=c(endog_to_rem,hks[!hks%in%hksToKeep])
keep=!ns.anno$Name%in%exclude
ns.data<-ns.data[keep,]
dim(ns.data)
ns.anno<-ns.anno[keep,]
dim(ns.anno)


## Add housekeeping normalization

nsNorm.cc.bg.hk<-NanoStringNorm(ns.data,
                                anno = ns.anno,
                                header = ns.header,
                                CodeCount = "geo.mean",
                                Background = "mean",
                                SampleContent = "housekeeping.sum",
                                round.values = T,take.log = T,traits = traits)

## Density plots

Plot.NanoStringNorm(nsNorm.cc.bg,plot.type = 'cv',title = F)
title("Normalized for Just Neg/Pos Controls")
Plot.NanoStringNorm(nsNorm.cc.bg.hk,plot.type = 'cv',title = F)
title("Adding Ref. Gene Norm")


## Check expression of replicates accross cartridges - I included replicates on each of the cartridges. i included 3, but one had quite low expression so I used the mean of 2

geomean<-function(x){exp(mean(log(x),na.rm=T))}
reps<-data.frame(cart1=apply(nsNorm.cc.bg.hk$normalized.data[,c("s7681a","s7667a")],1,geomean),
                 cart2=apply(nsNorm.cc.bg.hk$normalized.data[,c("s7681b","s7667b")],1,geomean),
                 cart3=apply(nsNorm.cc.bg.hk$normalized.data[,c("s7681c","s7667c")],1,geomean),
                 cart4=apply(nsNorm.cc.bg.hk$normalized.data[,c("s7681d","s7667d")],1,geomean),
                 row.names = rownames(nsNorm.cc.bg.hk$normalized.data))


boxplot(reps)
title("Boxplots of expression in technical replicates accross cartridges")

### Use mean of the replicates
RepMean<-function(x){
  # exclude Snca (human) which is 0 in WT
  x[!rownames(x)%in%"Snca","s7653a"]<-
    rowMeans(x[!rownames(x)%in%"Snca",c("s7653a","s7653b","s7653c","s7653d")], na.rm=T)
  x[!rownames(x)%in%"Snca","s7681a"]<-
    rowMeans(x[!rownames(x)%in%"Snca",c("s7681a","s7681b","s7681c","s7681d")], na.rm=T)
  x[!rownames(x)%in%"Snca","s7667a"]<-
    rowMeans(x[!rownames(x)%in%"Snca",c("s7667a","s7667b","s7667c","s7667d")],na.rm=T)
  x<-x[,!colnames(x)%in%c(
    "s7653b","s7653c","s7653d","s7681b","s7681c","s7681d","s7667b","s7667c","s7667d")]
  colnames(x)[10:12]<-c("s7653","s7681","s7667")
  return(x)
}

nsNorm<-nsNorm.cc.bg.hk$normalized.data[,-c(1:3)]
nsNorm<-RepMean(nsNorm)
samps<-Nano_SampleInfo[!Nano_SampleInfo$SampleNames%in%c("s7653b","s7653c","s7653d","s7681b",
                                                         "s7681c","s7681d","s7667b","s7667c","s7667d"),]
samps$SampleNames[10:12]<-c("s7653","s7681","s7667")
rownames(samps)<-samps$SampleNames
nsNorm<-nsNorm[rownames(nsNorm)%in%genes,]


save(nsNorm,samps,genes,hks,endog, file = "NSdata.rda")

# correlations with RNAseq--------
names(nsNorm)<-paste(gsub("s","",names(nsNorm)),"dc",sep = "_")

nsNorm2<-nsNorm[!rownames(nsNorm)%in%"Snca",]

eset<-readRDS("../../1-RNAseqWorkflow/C_QCandNormalization/data/rawCountEset.rds")

info<-fData(eset)
nsinfo<-info[info$Gene.name%in%rownames(nsNorm2),]
nsinfo<-nsinfo[order(match(nsinfo$Gene.name,rownames(nsNorm2))),]

rownames(nsNorm2)<-rownames(nsinfo)

sortByVar<-function(x){
  x<-x[order(apply(x,1,var),decreasing = T),]
  return(x)
}

sortByMean<-function(x){
  x<-x[order(apply(x,1,mean),decreasing = T),]
  return(x)
}

keepRows.value.in.n<-function(x,value,n=NULL,percent=0.5){
  #percent as decimel
  temp<-as.data.frame(x)
  temp$flag<-0
  for (i in 1:ncol(x)){
    temp$flag<-temp$flag+(temp[,i]>=value)
  }
  N<-ifelse(is.null(n),percent*ncol(x),n)
  keep=temp$flag>=N
  return(x[keep,])
}

ns2<-keepRows.value.in.n(nsNorm2,value = log2(50),percent = .3) #53 rows
topByMean<-rownames(sortByMean(ns2))
topByVar<-rownames(sortByVar(ns2))
keep<-union(topByMean[1:20],topByVar[1:20])

ns3<-ns2

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

library(Hmisc)
library(Biobase)
qsdat<-exprs(readRDS("../../1-RNAseqWorkflow/C_QCandNormalization/dc/data/CorrectedEset.rds"))

rows<-rownames(qsdat)[rownames(qsdat)%in%rownames(ns3)]
ns<-as.data.frame(t(ns3[rows,colnames(qsdat)]))
qs<-as.data.frame(t(qsdat[rows,]))

cors<-rcorr(as.matrix(qs),as.matrix(ns))
cors<-flattenCorrMatrix(cors$r, cors$P)
cors<-cors[cors$row==cors$column,]
cors$Symb<-nsinfo[cors$row,"Gene.name"]
cors$padj<-p.adjust(cors$p,method = "BH")
rownames(cors)<-cors$row
cors$avgCount<-rowMeans(ns3[rownames(cors),])
cors$var<-apply(ns3[rownames(cors),],1,var)

corsTopVar<-cors[cors$var>quantile(cors$var)[4],]

cors2<-cors[order(cors$p),c("Symb","cor","p","padj")]
write.csv(cors2,"nsCors.csv",row.names = F)
names(ns)<-nsinfo[rownames(ns3[rows,colnames(qsdat)]),"Gene.name"]
names(qs)<-names(ns)
library(tidyverse)
ns.l<-pivot_longer(rownames_to_column(ns,"Sample"),cols = c(2:ncol(ns)),names_to = "Gene",values_to = "nCounter")
qs.l<-pivot_longer(rownames_to_column(qs,"Sample"),cols = c(2:ncol(ns)),names_to = "Gene",values_to = "QuantSeq")
nsqs<-ns.l
nsqs$QuantSeq<-qs.l$QuantSeq


ggplot(nsqs[nsqs$Gene%in%cors2$Symb[cors2$padj<0.1],], aes(nCounter,QuantSeq,color=Gene))+
  geom_point()+
  geom_smooth(method = "lm",show.legend = F)+
  theme_bw()
ggsave("nsCorplot.pdf",width = 6.5,height = 4)

# transgene expression-----------

x=t(nsNorm)
pheno.grp<-samps
pheno.grp$time<-paste0("t",pheno.grp$t,"m")
nsTab<-cbind(Sample=colnames(nsNorm),
             gene=rep("SNCA",nrow(x)),
             Expr=x[,"Snca"],
             tissue=rep("Colon",nrow(x)),
             grp=paste0("dc",pheno.grp$t,"m"),
             model=rep("nCounter",nrow(x)),
             GT=as.character(pheno.grp$GT),
             Time=pheno.grp$time,
             SubGrp=rep("All",nrow(x))
)

nsTab<-as.data.frame(nsTab)
nsTab$Expr<-as.numeric(nsTab$Expr)
nsTab$Genotype<-gsub("Hem","Thy1-haSyn",nsTab$GT)
nsTab$Age[nsTab$Time=="t1m"]<-"1 month"
nsTab$Age[nsTab$Time=="t3m"]<-"3 months"

t.test(nsTab$Expr[nsTab$Genotype=="ASO"],
       nsTab$Expr[nsTab$Genotype=="WT"])
# t = 8.2826, df = 36.386, p-value = 6.802e-10
mean(nsTab$Expr[nsTab$Genotype=="ASO"])/mean(nsTab$Expr[nsTab$Genotype=="WT"])

t.test(nsTab$Expr[nsTab$Genotype=="ASO"&nsTab$Time=="t1m"],
       nsTab$Expr[nsTab$Genotype=="WT"&nsTab$Time=="t1m"])
# t = 6.889, df = 16.126, p-value = 3.492e-06
mean(nsTab$Expr[nsTab$Genotype=="ASO"&nsTab$Time=="t1m"])/mean(nsTab$Expr[nsTab$Genotype=="WT"&nsTab$Time=="t1m"])

# FC=8.2

t.test(nsTab$Expr[nsTab$Genotype=="ASO"&nsTab$Time=="t3m"],
       nsTab$Expr[nsTab$Genotype=="WT"&nsTab$Time=="t3m"])
#t = 4.4288, df = 11.842, p-value = 0.0008492
mean(nsTab$Expr[nsTab$Genotype=="ASO"&nsTab$Time=="t3m"])/mean(nsTab$Expr[nsTab$Genotype=="WT"&nsTab$Time=="t3m"])
# FC=3.4
library(ggplot2)
ggplot(nsTab[nsTab$gene=="SNCA",],
       aes(Age,Expr,color=Genotype))+
  geom_boxplot(outlier.size = 0,position = position_dodge(width = 1))+
  geom_dotplot(binaxis = "y",stackdir = "center",position = position_dodge(width = 1),binwidth = .1,aes(fill=Genotype))+
  theme_bw()+
  scale_fill_manual(values = c("#2774AE","#FFD100"))+
  scale_color_manual(values = c("#003B5C","#FFB81C"))+
  ylab("colon hSNCA (norm. expr.)")+
  xlab("Timepoint")+
  annotate("text",x = 1,y=7,label="p<0.005",hjust="left")+
  annotate("text",x = 2,y=7,label="p<0.005",hjust="left")


