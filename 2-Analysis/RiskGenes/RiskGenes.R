

# n.b. issue with biomaRt is preventing full reproducibility from raw data. The work-around which sets host = "https://dec2021.archive.ensembl.org/ results in one gene that is mapped from human to mouse differently (ENSG00000255439, VKORC1 And PRSS53 Readthrough is not mapped to mouse Vkorc1). The original annotation files can be loaded from /data.


library(tidyverse)
rm(list=ls())
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))


# 1. setup----------------
# create the annotation file to map human to mouse or load data used for analysis below in step 2

gwasStudies<-read.delim("data/gwasStudies.tsv",sep = "\t")

gwasAssoc<-read.delim("data/gwasAssociations.tsv",sep = "\t")

# setup PD-------
pdkeep<-grepl("Parkinson",gwasAssoc$MAPPED_TRAIT,ignore.case = T)
pdSnps<-gwasAssoc[pdkeep,]

pdSnps.intergen<-pdSnps[pdSnps$UPSTREAM_GENE_ID!="",]
pdSnps.intergen.l<-pivot_longer(pdSnps.intergen,
                                cols = c("UPSTREAM_GENE_ID","DOWNSTREAM_GENE_ID"),
                                names_to = "Type",values_to = "ENSG")
pdSnps.intergen.l$Type[pdSnps.intergen.l$Type=="UPSTREAM_GENE_ID"]<-"Upstream"
pdSnps.intergen.l$Type[pdSnps.intergen.l$Type=="DOWNSTREAM_GENE_ID"]<-"Downstream"

pdSnps.gene<-pdSnps[pdSnps$UPSTREAM_GENE_ID=="",]
snpgenes0<-strsplit(pdSnps.gene$SNP_GENE_IDS,", ")
max(unlist(lapply(snpgenes0,length))) #7
newnames<-paste("ID",seq(1,7,1),sep = ".")
snpgenes<-separate(pdSnps.gene,col = "SNP_GENE_IDS",sep = ", ",into =newnames,fill = "right" )

pdSnps.gene.l<-pivot_longer(snpgenes,
                            cols = all_of(newnames),
                            names_to = "Type",values_to = "ENSG")

pdSnps.gene.l<-pdSnps.gene.l[!is.na(pdSnps.gene.l$ENSG),]
pdSnps.gene.l$Type<-"Snp"

pdSnps.gene.l<-pdSnps.gene.l[,intersect(names(pdSnps.gene.l),names(pdSnps.intergen.l))]
pdSnps.intergen.l<-pdSnps.intergen.l[,intersect(names(pdSnps.gene.l),names(pdSnps.intergen.l))]

pdsnps.l<-rbind(pdSnps.gene.l,pdSnps.intergen.l[,names(pdSnps.gene.l)])


library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")

humGenes <- unique(pdsnps.l$ENSG)


# map human to mouse 

annotPD = getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"), 
               filters = "ensembl_gene_id", values = humGenes , 
               mart = human, attributesL = c("ensembl_gene_id","mgi_symbol"), 
               martL = mouse,uniqueRows = T)


annotPD[annotPD$Gene.stable.ID.1%in%annotPD$Gene.stable.ID.1[duplicated(annotPD$Gene.stable.ID.1)],]

pdsnps.l$ENSG2<-pdsnps.l$ENSG
pdsnps.l$ENSG2[pdsnps.l$ENSG=="ENSG00000288520"]<-"ENSG00000153721"
pdsnps.l$ENSG2[pdsnps.l$ENSG=="ENSG00000272414"]<-"ENSG00000189157"

annotPD<-annotPD[annotPD$Gene.stable.ID%in%pdsnps.l$ENSG2,]

pdsnp2<-NULL
for(mID in annotPD$Gene.stable.ID.1){
  hID=annotPD$Gene.stable.ID[annotPD$Gene.stable.ID.1==mID]
  
  snpRows<-pdsnps.l[pdsnps.l$ENSG2==hID,]
  snpRows$ENSGmus<-mID
  pdsnp2<-rbind(pdsnp2,snpRows)
}

saveRDS(pdsnp2,file = "data/PDsnpData.rds")
saveRDS(annotPD,file = "data/annotPD.rds")
# write.csv(pdsnp2,"PDsnpData.csv",row.names = F,na="")

PDkeyByMouseID<-data.frame(ENSG=unique(pdsnp2$ENSGmus),
                         PD.PMIDs=character(nrow(annotPD)))
for(i in 1:nrow(PDkeyByMouseID)){
  PDkeyByMouseID$PD.PMIDs[i]<-paste(pdsnp2$PUBMEDID[pdsnp2$ENSGmus==PDkeyByMouseID$ENSG[i]],collapse = ";")
  
}
saveRDS(PDkeyByMouseID,file = "data/PDkeyByMouseID.rds")



# setup IBD--------
ibdkeep<-grepl("Inflammatory Bowel Disease",gwasAssoc$MAPPED_TRAIT,ignore.case = T)|
  grepl("Ulcerative Colitis",gwasAssoc$MAPPED_TRAIT,ignore.case = T)|
  grepl("Crohn's",gwasAssoc$MAPPED_TRAIT,ignore.case = T)
IBDSnps<-gwasAssoc[ibdkeep,]

IBDSnps.intergen<-IBDSnps[IBDSnps$UPSTREAM_GENE_ID!="",]
IBDSnps.intergen.l<-pivot_longer(IBDSnps.intergen,
                                cols = c("UPSTREAM_GENE_ID","DOWNSTREAM_GENE_ID"),
                                names_to = "Type",values_to = "ENSG")
IBDSnps.intergen.l$Type[IBDSnps.intergen.l$Type=="UPSTREAM_GENE_ID"]<-"Upstream"
IBDSnps.intergen.l$Type[IBDSnps.intergen.l$Type=="DOWNSTREAM_GENE_ID"]<-"Downstream"

IBDSnps.gene<-IBDSnps[IBDSnps$UPSTREAM_GENE_ID=="",]
snpgenes0<-strsplit(IBDSnps.gene$SNP_GENE_IDS,", ")
max(unlist(lapply(snpgenes0,length))) #14
newnames<-paste("ID",seq(1,14,1),sep = ".")
snpgenes<-separate(IBDSnps.gene,col = "SNP_GENE_IDS",sep = ", ",into =newnames,fill = "right" )

IBDSnps.gene.l<-pivot_longer(snpgenes,
                            cols = all_of(newnames),
                            names_to = "Type",values_to = "ENSG")

IBDSnps.gene.l<-IBDSnps.gene.l[!is.na(IBDSnps.gene.l$ENSG),]
IBDSnps.gene.l$Type<-"Snp"

IBDSnps.gene.l<-IBDSnps.gene.l[,intersect(names(IBDSnps.gene.l),names(IBDSnps.intergen.l))]
IBDSnps.intergen.l<-IBDSnps.intergen.l[,intersect(names(IBDSnps.gene.l),names(IBDSnps.intergen.l))]

ibdsnps.l<-rbind(IBDSnps.gene.l,IBDSnps.intergen.l[,names(IBDSnps.gene.l)])

# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

humGenes <- unique(ibdsnps.l$ENSG)


annot = getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"), 
               filters = "ensembl_gene_id", values = humGenes , 
               mart = human, attributesL = c("ensembl_gene_id","mgi_symbol"), 
               martL = mouse,uniqueRows = T)

annotIBD<-annot
# make mapping unique - no duplicated mouse IDs
annotIBD[annotIBD$Gene.stable.ID.1%in%annotIBD$Gene.stable.ID.1[duplicated(annotIBD$Gene.stable.ID.1)],]
ibdsnps.l$ENSG2<-ibdsnps.l$ENSG
ibdsnps.l$ENSG2[ibdsnps.l$ENSG=="ENSG00000283321"]<-"ENSG00000106546" #Ahr
ibdsnps.l$ENSG2[ibdsnps.l$ENSG=="ENSG00000272305"]<-"ENSG00000163933" #Rft1
ibdsnps.l$ENSG2[ibdsnps.l$ENSG=="ENSG00000244255"]<-"ENSG00000243649" #Cfb
ibdsnps.l$ENSG2[ibdsnps.l$ENSG=="ENSG00000026036"]<-"ENSG00000258366" #Rtel1
ibdsnps.l$ENSG2[ibdsnps.l$ENSG=="ENSG00000281883"]<-"ENSG00000179630" #Lacc1
ibdsnps.l$ENSG2[ibdsnps.l$ENSG=="ENSG00000267303"]<-"ENSG00000161847" #Raver1
ibdsnps.l$ENSG2[ibdsnps.l$ENSG=="ENSG00000273154"]<-"ENSG00000197114" #Zgpat
ibdsnps.l$ENSG2[ibdsnps.l$ENSG=="ENSG00000285708"]<-"ENSG00000114861" #Foxp1
ibdsnps.l$ENSG2[ibdsnps.l$ENSG=="ENSG00000248672"]<-"ENSG00000054219" #Ly75



annotIBD<-annotIBD[annotIBD$Gene.stable.ID%in%ibdsnps.l$ENSG2,]

ibdsnp2<-NULL
for(mID in annotIBD$Gene.stable.ID.1){
  hID=annotIBD$Gene.stable.ID[annotIBD$Gene.stable.ID.1==mID]
  
  snpRows<-ibdsnps.l[ibdsnps.l$ENSG2==hID,]
  snpRows$ENSGmus<-mID
  ibdsnp2<-rbind(ibdsnp2,snpRows)
}

saveRDS(ibdsnp2,file = "data/IBDsnpData.rds")

# write.csv(ibdsnp2,"IBDsnpData.csv",row.names = F,na="")

IBDkeyByMouseID<-data.frame(ENSG=unique(ibdsnp2$ENSGmus),
                         IBD.PMIDs=character(nrow(annotIBD)))
for(i in 1:nrow(IBDkeyByMouseID)){
  IBDkeyByMouseID$IBD.PMIDs[i]<-paste(ibdsnp2$PUBMEDID[ibdsnp2$ENSGmus==IBDkeyByMouseID$ENSG[i]],collapse = ";")
  
}
saveRDS(IBDkeyByMouseID, file = "data/IBDkeyByMouseID.rds")
allSnpGenes<-union(IBDkeyByMouseID$ENSG,PDkeyByMouseID$ENSG)

# 2. load data-----------

gi.dc<-readRDS("../WGCNA/Colon/geneINFOall.rds") # gene info from WGCNA
gi.str<-readRDS("../WGCNA/Striatum/geneINFOallwithGrey.rds") # gene info from WGCNA
PDkeyByMouseID<-readRDS("data/PDkeyByMouseID-original.rds")
IBDkeyByMouseID<-readRDS("data/IBDkeyByMouseID-original.rds")
ibdsnpdat<-readRDS("data/IBDsnpData-original.rds")
pdsnpdat<-readRDS("data/PDsnpData-original.rds")



# 3. merge with gene info--------------
allSnpGenes<-union(IBDkeyByMouseID$ENSG,PDkeyByMouseID$ENSG)

gi.dc<-gi.dc[,c(1:2,4:6)]
gi.str<-gi.str[,c(1:2,4:6)]

gi.dc<-gi.dc[!gi.dc$module%in%"grey",]
gi.str<-gi.str[!gi.str$module%in%"grey",]

names(gi.dc)<-c( "ENSG","geneSymbol","Name","module.dc","MM.dc")

names(gi.str)<-c( "ENSG","geneSymbol","Name","module.str","MM.str")



gi.dc.names<-gi.dc[,c(1:3)]
gi.str.names<-gi.str[,c(1:3)]
rownames(gi.dc.names)<-NULL
rownames(gi.str.names)<-NULL

allnames<-rbind(gi.dc.names,gi.str.names)
allnames<-allnames[!duplicated(allnames$ENSG),]
snpNames<-allnames[allnames$ENSG%in%allSnpGenes,]

snpTable<-left_join(snpNames,PDkeyByMouseID,"ENSG")
snpTable<-left_join(snpTable,IBDkeyByMouseID,"ENSG")

rownames(gi.dc)<-NULL
rownames(gi.str)<-NULL
snpTable<-left_join(snpTable,gi.dc[,c(1,4:5)],"ENSG")
snpTable<-left_join(snpTable,gi.str[,c(1,4:5)],"ENSG")


pdsnpdat2<-pdsnpdat[pdsnpdat$ENSGmus%in%snpTable$ENSG,]
ibdsnpdat2<-ibdsnpdat[ibdsnpdat$ENSGmus%in%snpTable$ENSG,]


# 3. clean up ------------------

# abbreviated
unique(ibdsnpdat$DISEASE.TRAIT)

ibd.eq<-c("Chronic inflammatory diseases (ankylosing spondylitis, Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)",
          "Inflammatory bowel disease (early onset)","Pediatric autoimmune diseases",
          "Ulcerative colitis or Crohn's disease" 
)

cd.eq<-c("Crohn's disease and psoriasis" ,"Crohn's disease (need for surgery)" ,
         "Crohn's disease (indolent vs progressive)",
         "Crohn's disease (time to progression)" , 
         "Crohn's disease (time to first abdominal surgery)",
         "Poor prognosis in Crohn's disease"
)

exclude<-c("Thiopurine-induced leukopenia in inflammatory bowel disease" ,
           "Thiopurine-induced pancreatitis in inflammatory bowel disease" ,
           "Bone mineral density (lumbar spine) in inflammatory bowel disease",
           "Bone mineral density (femoral neck) in inflammatory bowel disease",
           "Thiopurine-induced alopecia in inflammatory bowel disease",
           "Thiopurine-induced myelosuppression in inflammatory bowel disease",
           "Paneth cell defects in Crohn's disease",
           "Pyoderma gangrenosum in inflammatory bowel disease","Crohn's disease and celiac disease",
           "Crohn's disease-related phenotypes",
           "Thiopurine-induced digestive symptoms in inflammatory bowel disease",
           "Infliximab-induced mucosal healing in Crohn's disease",
           "Erythema nodosum in inflammatory bowel disease",
           "Thiopurine-induced severe leukopenia in inflammatory bowel disease",
           "Thiopurine-induced acute severe leukopenia in inflammatory bowel disease",
           "Thiopurine-induced severe alopecia in inflammatory bowel disease",
           "Thiopurine-induced leukopenia in inflammatory bowel disease (conditioned on rs116855232)"
)
ibdsnpdat2<-ibdsnpdat[!ibdsnpdat$DISEASE.TRAIT%in%exclude,]
ibdsnpdat2$trait2<-ibdsnpdat2$DISEASE.TRAIT
ibdsnpdat2$trait2[ibdsnpdat2$trait2%in%ibd.eq]<-"IBD"
ibdsnpdat2$trait2[ibdsnpdat2$trait2%in%cd.eq]<-"CD"

ibdsnpdat2$trait2<-gsub("inflammatory bowel disease","IBD",ibdsnpdat2$trait2,ignore.case = T)
ibdsnpdat2$trait2<-gsub("Crohn's disease","CD",ibdsnpdat2$trait2,ignore.case = T)
ibdsnpdat2$trait2<-gsub("Ulcerative colitis","UC",ibdsnpdat2$trait2,ignore.case = T)

unique(ibdsnpdat2$trait2)

unique(pdsnpdat$DISEASE.TRAIT)
exclude<-c("Parkinsonism in frontotemporal lobe dementia")

pdsnpdat2<-pdsnpdat[!pdsnpdat$DISEASE.TRAIT%in%exclude,]
pdsnpdat2$trait2<-pdsnpdat2$DISEASE.TRAIT
pdsnpdat2$trait2[pdsnpdat2$trait2=="Parkinson's disease in GBA mutation carriers"]<-"PD (GBA)"
pdsnpdat2$trait2[pdsnpdat2$trait2=="Parkinson's disease (motor and cognition)"]<-"PD"

pdsnpdat2$trait2[pdsnpdat2$trait2=="Parkinson disease and lewy body pathology"]<-"PD/LBD"
pdsnpdat2$trait2[pdsnpdat2$trait2=="Parkinson's disease or first degree relation to individual with Parkinson's disease"]<-"PD or fam hx PD"

pdsnpdat2$trait2<-gsub("Parkinson's disease","PD",pdsnpdat2$trait2,ignore.case = T)
unique(pdsnpdat2$trait2)

saveRDS(ibdsnpdat2, file = "IBDsnpData-refined.rds")
saveRDS(pdsnpdat2, file = "PDsnpData-refined.rds")





ibdsnpsbytrait<-split(ibdsnpdat2,ibdsnpdat2$DISEASE.TRAIT)
pdsnpsbytrait<-split(pdsnpdat2,pdsnpdat2$DISEASE.TRAIT)

for(i in 1:nrow(snpTable)){
  if(snpTable$ENSG[i]%in%ibdsnpdat2$ENSGmus){
    cit<-NULL
    
    for(trait in names(ibdsnpsbytrait)){
      trait.df<-ibdsnpsbytrait[[trait]]
      trait.cits<-unique(trait.df$PUBMEDID[trait.df$ENSGmus==snpTable$ENSG[i]])
      
      if(length(trait.cits)>0){
        cit0<-paste0(trait," [", 
                     paste(trait.cits,collapse = ";") ,
                     "]")
        cit<-c(cit,cit0)
      }
      
    }
    snpTable$IBD.cits[i]<-paste(cit,collapse = ", ")
    
  }else{snpTable$IBD.cits[i]<-""}
  
}


for(i in 1:nrow(snpTable)){
  if(snpTable$ENSG[i]%in%pdsnpdat2$ENSGmus){
    cit<-NULL
    
    for(trait in names(pdsnpsbytrait)){
      trait.df<-pdsnpsbytrait[[trait]]
      trait.cits<-unique(trait.df$PUBMEDID[trait.df$ENSGmus==snpTable$ENSG[i]])
      
      if(length(trait.cits)>0){
        cit0<-paste0(trait," [", 
                     paste(trait.cits,collapse = ";") ,
                     "]")
        cit<-c(cit,cit0)
      }
      
    }
    snpTable$PD.cits[i]<-paste(cit,collapse = ", ")
    
  }else{snpTable$PD.cits[i]<-""}
  
}


saveRDS(snpTable,file = "data/snpTable.rds")


# 4. stats ---------------

snpTable2<-snpTable
snpTable2[snpTable2==""]<-NA
snpTable2$IBDsnp<-ifelse(is.na(snpTable2$IBD.cits),"No","Yes")
snpTable2$PDsnp<-ifelse(is.na(snpTable2$PD.cits),"No","Yes")

snpTable2.dc<-snpTable2[,c("ENSG","geneSymbol","Name" , 
                           "module.dc", "MM.dc" ,
                           "module.str","MM.str",
                           "IBDsnp","PDsnp")]
names(snpTable2.dc)[8:9]<-c("IBD","PD")
snpTable2.dc$module<-snpTable2.dc$module.dc
snpTable2.dc$MM<-snpTable2.dc$MM.dc

snpTable2.dc.alldc<-snpTable2.dc[!is.na(snpTable2.dc$MM.dc),]

snpTable2.str<-snpTable2[,c("ENSG","geneSymbol","Name" , 
                            "module.dc", "MM.dc" ,
                            "module.str","MM.str",
                            "IBDsnp","PDsnp")]
names(snpTable2.str)[8:9]<-c("IBD","PD")
snpTable2.str$module<-snpTable2.str$module.str
snpTable2.str$MM<-snpTable2.str$MM.str

snpTable2.str.allstr<-snpTable2.str[!is.na(snpTable2.str$MM.str),]

snptablists<-list(dc=snpTable2.dc.alldc,
                  str=snpTable2.str.allstr)


modlists<-list(dc=union(c("violet","blue","brown","darkmagenta",
                          "yellow","green","cyan","darkgreen",
                          "darkgrey","lightcyan","sienna3",
                          "greenyellow","pink"),
                        c("violet","brown","darkmagenta",
                          "yellow","green","cyan","darkgreen","sienna3",
                          "greenyellow","pink")),
               str=union(c("black", "cyan", "magenta", "red"),
                         c("cyan", "magenta", "red")))


gilist<-list(dc=readRDS("../WGCNA/Colon/geneINFOall.rds"),
             str=readRDS("../WGCNA/Striatum/geneINFOallwithGrey.rds"))

RiskGenes=list(PD=readRDS("PDsnpData-refined.rds")$ENSGmus,
               IBD=readRDS("IBDsnpData-refined.rds")$ENSGmus)

allgenes.list=list(dc=gilist$dc$ENSG,
                   str=gilist$str$ENSG)




snprestabByModule=NULL
snprestab=NULL
for(Disease in c("PD","IBD")){
  for(Tissue in c("dc","str")){
    keymods=modlists[[Tissue]]
    setgenes=RiskGenes[[Disease]]
    allgenes=allgenes.list[[Tissue]]
    setgenesinallgenes<-intersect(setgenes,allgenes)
    n.geneset<-length(setgenesinallgenes)
    n.all<-length(allgenes)
    gitab<-gilist[[Tissue]]
    gitab<-gitab[gitab$module!="grey",]
    gitab<-gitab[gitab$ENSG%in%allgenes,]
    gitab$ASOmodule<-ifelse(gitab$module%in%keymods,"Yes","No")
    gitab[,Disease]<-ifelse(gitab$ENSG%in%setgenes,"Yes","No")
    for(MM in c(0.2,0.5)){
      gitab<-gitab[gitab$MM>=MM,]
      hyper<-NULL
      for(mod in unique(gitab$module)){
        modgenes<-gitab$ENSG[gitab$moduleColor==mod]
        n.module<-length(modgenes)
        n.geneset<-length(setgenesinallgenes)
        n.all<-length(allgenes)
        overlap<-intersect(setgenesinallgenes,modgenes)
        n.overlap<-length(overlap)
        pval<-phyper(n.overlap-1, n.module, n.all-n.module, 
                     n.geneset,lower.tail= FALSE)
        
        overlapSymbs<-paste(gitab$geneSymbol[gitab$ENSG%in%overlap],collapse = "/")
        hyper<-rbind(hyper,
                     c(Tissue,Disease,
                       mod,n.overlap,pval,
                       n.geneset,n.module,overlapSymbs,paste0("MM-",MM)))
      }
      
      hyper<-as.data.frame(hyper)
      names(hyper)<-c("Tissue","Disease",
                      "Module","Overlap","P",
                      "SetSize","ModSize",
                      "CommonGenes","MMthresh")
      hyper$ASOmodule<-ifelse(hyper$Module%in%keymods,"Yes","No")
      f.genes<-fisher.test(table(gitab$ASOmodule,
                                 gitab[,Disease]))
      
      hyper$p.ft.genes<-f.genes$p.value
      f.modules<-fisher.test(table(hyper$ASOmodule,hyper$Overlap))
      hyper$p.ft.modules<-f.modules$p.value
      
      snprestabByModule<-rbind(snprestabByModule,hyper)
      snprestab0<-hyper[1,!names(hyper)%in%c("Module","Overlap","P",
                                             "SetSize",	"ModSize",	"CommonGenes",
                                             "ASOmodule")]
      hyper2<-hyper
      
      hyper2$Mod2<-paste0(hyper2$Module,
                          " (",
                          hyper2$SetSize,
                          ", ",
                          hyper2$ModSize,
                          ", ",
                          hyper2$Overlap,
                          ")")
      snprestab0$HyperMods<-paste(hyper2$Module[hyper2$P<0.05],collapse = "/")
      
      snprestab0$HyperModsSize<-paste(hyper2$Mod2[hyper2$P<0.05],collapse = "/")
      
      
      snprestab<-rbind(snprestab,snprestab0)
    }
    
  }
}


write.csv(snprestab,"snprestab.csv",row.names = F)
write.csv(snprestabByModule,"snprestabByModule.csv",row.names = F)


