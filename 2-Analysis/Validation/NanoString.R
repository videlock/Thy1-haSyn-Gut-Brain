
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))

# NanoStringNorm needs to be installed from source
install.packages(pkgs="NanoStringNorm_1.2.1.1.tar.gz", type="source", repos=NULL)

library(NanoStringNorm)


## load data


load("NSraw.rda")
Nano_SampleInfo<-readRDS("NanoSampleInfo.rds")



ns.data<-NanoString.mRNA$x[,-c(1:3)]
ns.anno <- NanoString.mRNA$x[,c(1:3)]
ns.header<-NanoString.mRNA$header
traits<-data.frame(gt=ifelse(Nano_SampleInfo$GT=="Hem",1,2),
                   time=ifelse(Nano_SampleInfo$Time=="4 wks",1,2))
rownames(traits)<-NULL

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
find_outlier <- function(x,mult=1.5) {
  return(x < quantile(x, .25,na.rm=T) - mult*IQR(x,na.rm=T) | x > quantile(x, .75,na.rm=T) + mult*IQR(x,na.rm=T))
}


ns<-NanoStringNorm(ns.data,
                    anno = ns.anno,
                   header = ns.header,
                   CodeCount = "geo.mean",
                   Background = "mean",
                   SampleContent = "housekeeping.sum",
                   round.values = T,take.log = F, traits = traits)

datNS<-ns$normalized.data[,-c(1:3)]
datNS<-RepMean(datNS)
samps<-Nano_SampleInfo[!Nano_SampleInfo$SampleNames%in%c("s7653b","s7653c","s7653d","s7681b",
                                                         "s7681c","s7681d","s7667b","s7667c","s7667d"),]
samps$SampleNames[10:12]<-c("s7653","s7681","s7667")
rownames(samps)<-samps$SampleNames
samps$t<-ifelse(samps$Time=="4 wks",1,3)
samps$GTtime<-paste(samps$GT,samps$t,sep = "_")
datNS<-datNS[rownames(datNS)%in%genes,]


x=t(datNS)
pheno.grp<-samps
pheno.grp$time<-paste0("t",pheno.grp$t,"m")



nsTab<-cbind(Sample=colnames(datNS),
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

nsTab <- nsTab %>% 
  group_by(Age,GT) %>% 
  mutate(Expr2=ifelse(find_outlier(Expr),NA,Expr)) %>% 
  ungroup()


t.test(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"],
       nsTab$Expr2[nsTab$Genotype=="WT"])

mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"],na.rm=T)


mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t1m"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t1m"],na.rm=T)

mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t3m"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t3m"],na.rm=T)

nstab.hs<-nsTab


nsTab<-cbind(Sample=colnames(datNS),
             gene=rep("Snca",nrow(x)),
             Expr=x[,"Snca Mm"],
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


nsTab <- nsTab %>% 
  group_by(Age,GT) %>% 
  mutate(Expr2=ifelse(find_outlier(Expr),NA,Expr)) %>% 
  ungroup()


t.test(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"],
       nsTab$Expr2[nsTab$Genotype=="WT"])



mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"],na.rm=T)

t.test(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t1m"],
       nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t1m"])


mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t1m"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t1m"],na.rm=T)



t.test(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t3m"],
       nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t3m"])



mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t3m"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t3m"],na.rm=T)

nstab.mm<-nsTab


# combined--------

nstab<-rbind(nstab.hs,nstab.mm)


library(ggpubr)
library(rstatix)

nstab$gene2 <- ifelse(nstab$gene=="Snca",
                      "Snca (murine)",
                      "SNCA (transgene)")
nstab$gene2 <- factor(nstab$gene2,levels = c("SNCA (transgene)","Snca (murine)"))
p <- ggplot(nstab,
            aes(Age,Expr2))+
  geom_boxplot(outlier.size = 0,position = position_dodge(width = 1),
               aes(Age,Expr2,color=Genotype))+
  theme_bw()+
  geom_point(aes(Age,Expr2,group=Genotype,color=Genotype),position = position_dodge(width = 1),show.legend = F)+
  scale_color_manual(values = c("#003B5C","#FFB81C"))+
  ylab("colon aSyn (norm. exprs.)")+
  xlab("Timepoint")+
  facet_wrap(~gene2)+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold"))

stat.test <- nstab %>%
  group_by(Age,gene2) %>%
  wilcox_test(Expr2 ~ GT, ref.group = "WT") %>%
  add_significance(cutpoints = c(0,0.005,0.05,1),
                   symbols = c("**","*",""))

stat.test <- stat.test %>% add_xy_position(x = "Age",dodge = 1)


p + stat_pvalue_manual(stat.test)




write_rds(nstab, file = "nstab.rds")






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






ns<-NanoStringNorm(ns.data,
                   anno = ns.anno,
                   header = ns.header,
                   CodeCount = "geo.mean",
                   Background = "mean",
                   SampleContent = "housekeeping.sum",
                   round.values = T,take.log = F, traits = traits)

datNS<-ns$normalized.data[,-c(1:3)]
datNS<-RepMean(datNS)
samps<-Nano_SampleInfo[!Nano_SampleInfo$SampleNames%in%c("s7653b","s7653c","s7653d","s7681b",
                                                         "s7681c","s7681d","s7667b","s7667c","s7667d"),]
samps$SampleNames[10:12]<-c("s7653","s7681","s7667")
rownames(samps)<-samps$SampleNames
samps$t<-ifelse(samps$Time=="4 wks",1,3)
samps$GTtime<-paste(samps$GT,samps$t,sep = "_")
datNS<-datNS[rownames(datNS)%in%genes,]


x=t(datNS)
pheno.grp<-samps
pheno.grp$time<-paste0("t",pheno.grp$t,"m")
nsTab<-cbind(Sample=colnames(datNS),
             gene=rep("Snca",nrow(x)),
             Expr=x[,"Snca Mm"],
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


nsTab <- nsTab %>% 
  group_by(Age,GT) %>% 
  mutate(Expr2=ifelse(find_outlier(Expr),NA,Expr)) %>% 
  ungroup()


t.test(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"],
       nsTab$Expr2[nsTab$Genotype=="WT"])
# t = 2.642, df = 17.034, p-value = 0.0171




mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"],na.rm=T)
# 84.17

t.test(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t1m"],
       nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t1m"])
# t = 1.3152, df = 8.4441, p-value = 0.223



mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t1m"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t1m"],na.rm=T)

# 9

t.test(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t3m"],
       nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t3m"])
#t = 2.8419, df = 8, p-value = 0.02175



mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t3m"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t3m"],na.rm=T)
# inf


ggplot(nsTab[nsTab$gene=="Snca",],
       aes(Age,Expr2,color=Genotype))+
  geom_boxplot(outlier.size = 0,position = position_dodge(width = 1))+
  geom_dotplot(binaxis = "y",stackdir = "center",position = position_dodge(width = 1),
               binwidth = .3,aes(fill=Genotype))+
  theme_bw()+
  scale_fill_manual(values = c("#2774AE","#FFD100"))+
  scale_color_manual(values = c("#003B5C","#FFB81C"))+
  ylab("colon Snca (gm pos cont, mean bg, hk sum)")+
  xlab("Timepoint")+
  annotate("text",x = 1,y=10,label="p=0.22",hjust="left")+
  annotate("text",x = 2,y=10,label="p=0.022",hjust="left")

ggsave("ncounterSnca-mm-pc-gm-bg-mn-hksum.pdf",width = 4,height = 3)
nstab.mm<-nsTab


# combined--------


ggplot(nstab.hs[nstab.hs$gene=="SNCA",],
       aes(Age,Expr2,color=Genotype))+
  geom_boxplot(outlier.size = 0,position = position_dodge(width = 1))+
  geom_dotplot(binaxis = "y",stackdir = "center",position = position_dodge(width = 1),
               binwidth = 4,aes(fill=Genotype))+
  theme_bw()+
  scale_fill_manual(values = c("#2774AE","#FFD100"))+
  scale_color_manual(values = c("#003B5C","#FFB81C"))+
  ylab("colon hSNCA (norm. exprs.)")+
  xlab("Timepoint")+
  annotate("text",x = 1,y=110,label="p<0.005",hjust="left")+
  annotate("text",x = 2,y=75,label="p<0.005",hjust="left")+
  ggtitle("A. hSNCA Transgene Expression")

ggsave("ncounterSNCA-revised.pdf",width = 4,height = 3)

ggplot(nstab.mm[nstab.mm$gene=="Snca",],
       aes(Age,Expr2,color=Genotype))+
  geom_boxplot(outlier.size = 0,position = position_dodge(width = 1))+
  geom_dotplot(binaxis = "y",stackdir = "center",position = position_dodge(width = 1),
               binwidth = .3,aes(fill=Genotype))+
  theme_bw()+
  scale_fill_manual(values = c("#2774AE","#FFD100"))+
  scale_color_manual(values = c("#003B5C","#FFB81C"))+
  ylab("colon Snca (norm. exprs.)")+
  xlab("Timepoint")+
  annotate("text",x = 1,y=10,label="p=0.22",hjust="left")+
  annotate("text",x = 2,y=10,label="p=0.022",hjust="left")

ggsave("ncounterSnca-mm-revised.pdf",width = 4,height = 3)

nstab<-rbind(nstab.hs,nstab.mm)
# nstab$outlier<-ifelse(nstab$Sample%in%c("s7745","s7644"),"outlier","")
ggplot(nstab,
       aes(Age,Expr2))+
  geom_boxplot(outlier.size = 0,position = position_dodge(width = 1),
               aes(Age,Expr2,color=Genotype))+
  # geom_dotplot(binaxis = "y",stackdir = "center",position = position_dodge(width = 1), aes(Age,Expr,fill=outlier),
  #              binwidth = 2)+
  theme_bw()+
  # geom_dotplot(data = nstab[nstab$Sample%in%c("s7745"),], binaxis = "y",stackdir = "center",
  #             binwidth = 2,color="green",aes(fill=Genotype), position = pos)+
  # scale_fill_manual(values = c("#2774AE","#FFD100"))+
  geom_point(aes(Age,Expr2,group=Genotype,color=Genotype),position = position_dodge(width = 1),show.legend = F)+
  scale_color_manual(values = c("#003B5C","#FFB81C"))+
  ylab("colon asyn (norm. exprs.)")+
  xlab("Timepoint")+
  facet_wrap(~gene)
ggsave("ncounterSNCASnca.pdf",width = 6,height = 4)


library(ggpubr)
library(rstatix)

nstab$gene2 <- ifelse(nstab$gene=="Snca",
                      "Snca (murine)",
                      "SNCA (transgene)")
nstab$gene2 <- factor(nstab$gene2,levels = c("SNCA (transgene)","Snca (murine)"))
p <- ggplot(nstab,
       aes(Age,Expr2))+
  geom_boxplot(outlier.size = 0,position = position_dodge(width = 1),
               aes(Age,Expr2,color=Genotype))+
  theme_bw()+
  geom_point(aes(Age,Expr2,group=Genotype,color=Genotype),position = position_dodge(width = 1),show.legend = F)+
  scale_color_manual(values = c("#003B5C","#FFB81C"))+
  ylab("colon aSyn (norm. exprs.)")+
  xlab("Timepoint")+
  facet_wrap(~gene2)+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold"))


write_rds(nstab, file = "nstab.rds")



p



stat.test <- nstab %>%
  group_by(Age,gene2) %>%
  wilcox_test(Expr2 ~ GT, ref.group = "WT") %>%
  add_significance(cutpoints = c(0,0.005,0.05,1),
                   symbols = c("**","*",""))

stat.test <- stat.test %>% add_xy_position(x = "Age",dodge = 1)


p + stat_pvalue_manual(stat.test)


ggsave("ncounterSNCASnca.pdf",width = 6,height = 4)













load("NSraw.rda")
Nano_SampleInfo<-readRDS("NanoSampleInfo.rds")
ns.data<-NanoString.mRNA$x[,-c(1:3)]
ns.anno <- NanoString.mRNA$x[,c(1:3)]
ns.header<-NanoString.mRNA$header
traits<-data.frame(gt=ifelse(Nano_SampleInfo$GT=="Hem",1,2),
                   time=ifelse(Nano_SampleInfo$Time=="4 wks",1,2))
# ns.anno$Name[ns.anno$Name=="Cacna2s3"]<-"Cacna2d3"
# keep=!colnames(ns.data)%in%"s7664"
# ns.data<-ns.data[,keep]
# ns.header<-ns.header[,keep]
# traits<-traits[keep,]
rownames(traits)<-NULL
# Nano_SampleInfo<-Nano_SampleInfo[keep,]
# 
# NanoString.mRNA$x<- cbind(ns.anno,ns.data)
# NanoString.mRNA$header <- ns.header
# save(NanoString.mRNA, endog, genes, hks, file = "NSraw.rda")
# saveRDS(Nano_SampleInfo, file = "NanoSampleInfo.rds")

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

# transgene expression-----------

find_outlier <- function(x,mult=1.5) {
  return(x < quantile(x, .25,na.rm=T) - mult*IQR(x,na.rm=T) | x > quantile(x, .75,na.rm=T) + mult*IQR(x,na.rm=T))
}



ns<-NanoStringNorm(ns.data,
                   anno = ns.anno,
                   header = ns.header,
                   CodeCount = "geo.mean",
                   Background = "mean",
                   SampleContent = "housekeeping.sum",
                   round.values = T,take.log = F, traits = traits)

datNS<-ns$normalized.data[,-c(1:3)]
datNS<-RepMean(datNS)
samps<-Nano_SampleInfo[!Nano_SampleInfo$SampleNames%in%c("s7653b","s7653c","s7653d","s7681b",
                                                         "s7681c","s7681d","s7667b","s7667c","s7667d"),]
samps$SampleNames[10:12]<-c("s7653","s7681","s7667")
rownames(samps)<-samps$SampleNames
samps$t<-ifelse(samps$Time=="4 wks",1,3)
samps$GTtime<-paste(samps$GT,samps$t,sep = "_")
datNS<-datNS[rownames(datNS)%in%genes,]


x=t(datNS)
pheno.grp<-samps
pheno.grp$time<-paste0("t",pheno.grp$t,"m")
nsTab<-cbind(Sample=colnames(datNS),
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

nsTab <- nsTab %>% 
  group_by(Age,GT) %>% 
  mutate(Expr2=ifelse(find_outlier(Expr),NA,Expr)) %>% 
  ungroup()


t.test(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"],
       nsTab$Expr2[nsTab$Genotype=="WT"])


mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"],na.rm=T)


mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t1m"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t1m"],na.rm=T)




mean(nsTab$Expr2[nsTab$Genotype=="Thy1-haSyn"&nsTab$Time=="t3m"],na.rm=T)/mean(nsTab$Expr2[nsTab$Genotype=="WT"&nsTab$Time=="t3m"],na.rm=T)


ggplot(nsTab[nsTab$gene=="SNCA",],
       aes(Age,Expr2,color=Genotype))+
  geom_boxplot(outlier.size = 0,position = position_dodge(width = 1))+
  geom_dotplot(binaxis = "y",stackdir = "center",position = position_dodge(width = 1),
               binwidth = 4,aes(fill=Genotype))+
  theme_bw()+
  scale_fill_manual(values = c("#2774AE","#FFD100"))+
  scale_color_manual(values = c("#003B5C","#FFB81C"))+
  ylab("colon hSNCA (gm pos cont, mean bg, hk sum)")+
  xlab("Timepoint")+
  annotate("text",x = 1,y=110,label="p<0.005",hjust="left")+
  annotate("text",x = 2,y=75,label="p<0.005",hjust="left")

ggsave("ncounterSNCA-revised.pdf",width = 4,height = 3)

nstab.hs<-nsTab






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


