
library(tidyverse)

# snps---------------
# BiocManager::install("gwascat")
# library(gwascat)

# setup----------------
gwasStudies<-read.delim("gwasStudies.tsv",sep = "\t")
# head(gwasStudies)

gwasAssoc<-read.delim("gwasAssociations.tsv",sep = "\t")
# head(gwasAssoc)

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
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

humGenes <- unique(pdsnps.l$ENSG)


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

saveRDS(pdsnp2,file = "PDsnpData.rds")
saveRDS(annotPD,file = "annotPD.rds")
write.csv(pdsnp2,"PDsnpData.csv",row.names = F,na="")

PDkeyByMouseID<-data.frame(ENSG=unique(pdsnp2$ENSGmus),
                         PD.PMIDs=character(nrow(annotPD)))
for(i in 1:nrow(PDkeyByMouseID)){
  PDkeyByMouseID$PD.PMIDs[i]<-paste(pdsnp2$PUBMEDID[pdsnp2$ENSGmus==PDkeyByMouseID$ENSG[i]],collapse = ";")
  
}
saveRDS(PDkeyByMouseID,file = "PDkeyByMouseID.rds")



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

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

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

saveRDS(ibdsnp2,file = "IBDsnpData.rds")

write.csv(ibdsnp2,"IBDsnpData.csv",row.names = F,na="")

IBDkeyByMouseID<-data.frame(ENSG=unique(ibdsnp2$ENSGmus),
                         IBD.PMIDs=character(nrow(annotIBD)))
for(i in 1:nrow(IBDkeyByMouseID)){
  IBDkeyByMouseID$IBD.PMIDs[i]<-paste(ibdsnp2$PUBMEDID[ibdsnp2$ENSGmus==IBDkeyByMouseID$ENSG[i]],collapse = ";")
  
}
saveRDS(IBDkeyByMouseID, file = "IBDkeyByMouseID.rds")
allSnpGenes<-union(IBDkeyByMouseID$ENSG,PDkeyByMouseID$ENSG)

# load -------
gi.dc<-readRDS("dc/geneINFOall.rds")
gi.str<-readRDS("str/geneINFOall.rds")
PDkeyByMouseID<-readRDS("PDkeyByMouseID.rds")
IBDkeyByMouseID<-readRDS("IBDkeyByMouseID.rds")
ibdsnpdat<-readRDS("IBDsnpData.rds")
pdsnpdat<-readRDS("PDsnpData.rds")

# merge with gene info--------------
allSnpGenes<-union(IBDkeyByMouseID$ENSG,PDkeyByMouseID$ENSG)

gi.dc<-gi.dc[,c(1:2,4:6)]
gi.str<-gi.str[,c(1:2,4:6)]

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




ibdPMIDs<-unique(ibdsnpdat$PUBMEDID[ibdsnpdat$ENSGmus%in%snpTable$ENSG])
PDPMIDs<-unique(pdsnpdat$PUBMEDID[pdsnpdat$ENSGmus%in%snpTable$ENSG])

allPMIDs<-union(ibdPMIDs,PDPMIDs)



# install.packages("RefManageR")
library(RefManageR)

SNPcitations<-GetPubMedByID(allPMIDs)
SNPcitTable<-as.data.frame(SNPcitations)
SNPcitTable2<-data.frame(PUBMEDID=SNPcitTable$eprint,citekey=paste0("@",rownames(SNPcitTable)))
SNPcitTable2$PUBMEDID<-as.integer(SNPcitTable2$PUBMEDID)

pdsnpdat2<-pdsnpdat[pdsnpdat$ENSGmus%in%snpTable$ENSG,]
pdsnpdat2<-left_join(pdsnpdat2,SNPcitTable2,"PUBMEDID")

ibdsnpdat2<-ibdsnpdat[ibdsnpdat$ENSGmus%in%snpTable$ENSG,]
ibdsnpdat2<-left_join(ibdsnpdat2,SNPcitTable2,"PUBMEDID")

for(i in 1:nrow(snpTable)){
  if(snpTable$ENSG[i]%in%pdsnpdat2$ENSGmus){
    snpTable$PDany.citation[i]<-paste0("[",
                                       paste(unique(pdsnpdat2$citekey[pdsnpdat2$ENSGmus==snpTable$ENSG[i]]),
                                             collapse = ";") ,
                                       "]"
                                       
    )
    
  }else{snpTable$PDany.citation[i]<-""}
  
}



for(i in 1:nrow(snpTable)){
  if(snpTable$ENSG[i]%in%ibdsnpdat2$ENSGmus){
    snpTable$IBDany.citation[i]<-paste0("[",
                                        paste(unique(ibdsnpdat2$citekey[ibdsnpdat2$ENSGmus==snpTable$ENSG[i]]),
                                              collapse = ";") ,
                                        "]"
                                        
    )
    
  }else{snpTable$IBDany.citation[i]<-""}
  
}

ibdsnpsbytrait<-split(ibdsnpdat2,ibdsnpdat2$DISEASE.TRAIT)
pdsnpsbytrait<-split(pdsnpdat2,pdsnpdat2$DISEASE.TRAIT)

for(i in 1:nrow(snpTable)){
  if(snpTable$ENSG[i]%in%ibdsnpdat2$ENSGmus){
    cit<-NULL
    
    for(trait in names(ibdsnpsbytrait)){
      trait.df<-ibdsnpsbytrait[[trait]]
      trait.cits<-unique(trait.df$citekey[trait.df$ENSGmus==snpTable$ENSG[i]])
      
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
      trait.cits<-unique(trait.df$citekey[trait.df$ENSGmus==snpTable$ENSG[i]])
      
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


# abrreviated
unique(ibdsnpdat2$DISEASE.TRAIT)

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
           "Erythema nodosum in inflammatory bowel disease"
)
ibdsnpdat3<-ibdsnpdat2[!ibdsnpdat2$DISEASE.TRAIT%in%exclude,]
ibdsnpdat3$trait2<-ibdsnpdat3$DISEASE.TRAIT
ibdsnpdat3$trait2[ibdsnpdat3$trait2%in%ibd.eq]<-"IBD"
ibdsnpdat3$trait2[ibdsnpdat3$trait2%in%cd.eq]<-"CD"

ibdsnpdat3$trait2<-gsub("inflammatory bowel disease","IBD",ibdsnpdat3$trait2,ignore.case = T)
ibdsnpdat3$trait2<-gsub("Crohn's disease","CD",ibdsnpdat3$trait2,ignore.case = T)
ibdsnpdat3$trait2<-gsub("Ulcerative colitis","UC",ibdsnpdat3$trait2,ignore.case = T)

unique(ibdsnpdat3$trait2)

unique(pdsnpdat2$DISEASE.TRAIT)
exclude<-c("Parkinsonism in frontotemporal lobe dementia")

pdsnpdat3<-pdsnpdat2[!pdsnpdat2$DISEASE.TRAIT%in%exclude,]
pdsnpdat3$trait2<-pdsnpdat3$DISEASE.TRAIT
pdsnpdat3$trait2[pdsnpdat3$trait2=="Parkinson's disease in GBA mutation carriers"]<-"PD (GBA)"
pdsnpdat3$trait2[pdsnpdat3$trait2=="Parkinson's disease (motor and cognition)"]<-"PD"

pdsnpdat3$trait2[pdsnpdat3$trait2=="Parkinson disease and lewy body pathology"]<-"PD/LBD"
pdsnpdat3$trait2[pdsnpdat3$trait2=="Parkinson's disease or first degree relation to individual with Parkinson's disease"]<-"PD or fam hx PD"

pdsnpdat3$trait2<-gsub("Parkinson's disease","PD",pdsnpdat3$trait2,ignore.case = T)
unique(pdsnpdat3$trait2)

ibdsnpsbytrait2<-split(ibdsnpdat3,ibdsnpdat3$trait2)
pdsnpsbytrait2<-split(pdsnpdat3,pdsnpdat3$trait2)

for(i in 1:nrow(snpTable)){
  if(snpTable$ENSG[i]%in%ibdsnpdat3$ENSGmus){
    cit<-NULL
    
    for(trait in names(ibdsnpsbytrait2)){
      trait.df<-ibdsnpsbytrait2[[trait]]
      trait.cits<-unique(trait.df$citekey[trait.df$ENSGmus==snpTable$ENSG[i]])
      
      if(length(trait.cits)>0){
        cit0<-paste0(trait," [", 
                     paste(trait.cits,collapse = ";") ,
                     "]")
        cit<-c(cit,cit0)
      }
      
    }
    snpTable$IBD.GWAS[i]<-paste(cit,collapse = ", ")
    
  }else{snpTable$IBD.GWAS[i]<-""}
  
}


for(i in 1:nrow(snpTable)){
  if(snpTable$ENSG[i]%in%pdsnpdat3$ENSGmus){
    cit<-NULL
    
    for(trait in names(pdsnpsbytrait2)){
      trait.df<-pdsnpsbytrait2[[trait]]
      trait.cits<-unique(trait.df$citekey[trait.df$ENSGmus==snpTable$ENSG[i]])
      
      if(length(trait.cits)>0){
        cit0<-paste0(trait," [", 
                     paste(trait.cits,collapse = ";") ,
                     "]")
        cit<-c(cit,cit0)
      }
      
    }
    snpTable$PD.GWAS[i]<-paste(cit,collapse = ", ")
    
  }else{snpTable$PD.GWAS[i]<-""}
  
}




saveRDS(snpTable,file = "snpTable.rds")