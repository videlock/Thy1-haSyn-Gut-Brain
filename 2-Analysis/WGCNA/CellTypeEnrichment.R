# uses marker lists downloaded from Panglao and CellMarker
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# fix the table------------
# mgenes<-getBM(attributes = c("mgi_symbol"),mart = mouse)
# h2genes<-mgenes$mgi_symbol[grep("H2-",mgenes$mgi_symbol)]
# paste(h2genes,collapse = ", ")
# 
# aldhhgenes<-mgenes$mgi_symbol[grep("Aldh",mgenes$mgi_symbol)]
# paste(aldhhgenes,collapse = ", ")
# 
# paste(mgenes$mgi_symbol[grep("Clec",mgenes$mgi_symbol)],collapse = ", ")
# 
# paste(mgenes$mgi_symbol[grep("S100",mgenes$mgi_symbol)],collapse = ", ")
# paste(mgenes$mgi_symbol[grep("Vegf",mgenes$mgi_symbol)],collapse = ", ")
# 
# paste(mgenes$mgi_symbol[grep("Ighd",mgenes$mgi_symbol)],collapse = ", ")
# 
# paste(mgenes$mgi_symbol[grep("Ighg",mgenes$mgi_symbol)],collapse = ", ")
# paste(mgenes$mgi_symbol[grep("Klra",mgenes$mgi_symbol)],collapse = ", ")
# paste(mgenes$mgi_symbol[grep("Tgfb",mgenes$mgi_symbol)],collapse = ", ")
# paste(mgenes$mgi_symbol[grep("Ly6",mgenes$mgi_symbol)],collapse = ", ")
# paste(mgenes$mgi_symbol[grep("Krt",mgenes$mgi_symbol)],collapse = ", ")
# paste(mgenes$mgi_symbol[grep("Ikzf",mgenes$mgi_symbol)],collapse = ", ")
# paste(mgenes$mgi_symbol[grep("Pitp",mgenes$mgi_symbol)],collapse = ", ")


# marker lists------
markersP<-read.delim("PanglaoDB_markers_27_Mar_2020.tsv")
markersP<-markersP[!markersP$species%in%"Hs",]
# unique(markersP$organ)

# write.table(markersP$official.gene.symbol,quote = F,sep = "\t",row.names = F,col.names = F,file = "panglaoHum.txt")

annot = getLDS(attributes = c("hgnc_symbol"), 
               filters = "hgnc_symbol", values = unique(markersP$official.gene.symbol), 
               mart = human, attributesL = c("mgi_symbol"), 
               martL = mouse, uniqueRows=T)



library(tidyverse)
names(annot)<-c("official.gene.symbol","mgi_symbol")
markersP<-left_join(markersP,annot,"official.gene.symbol")
markersP<-markersP[!is.na(markersP$mgi_symbol),]

markersP$CellType<-paste0(markersP$cell.type," (",markersP$organ,")")
markersP<-markersP[markersP$canonical.marker==1,]
markersP$db<-"Panglao"

P.tissues.dc.all<-c("GI tract","Blood","Connective tissue",
                    "Brain","Epithelium","Vasculature","Smooth muscle",
                    "Immune System",NA)
P.tissues.dc.gi<-c("GI tract")

P.tissues.str.all<-c("Brain","Blood","Immune System",
                     "Vasculature",
                     NA)
P.tissues.str.brain<-c("Brain")
markersP.dc<-markersP[markersP$organ%in%P.tissues.dc.all,]
markersP.dc.GIonly<-markersP[markersP$organ%in%P.tissues.dc.gi,]
markersP.str<-markersP[markersP$organ%in%P.tissues.str.all,]
markersP.str.BrainOnly<-markersP[markersP$organ%in%P.tissues.str.brain,]

markersListP.dc<-split(markersP.dc,markersP.dc$CellType)
markersListP.dc.GIonly<-split(markersP.dc.GIonly,markersP.dc.GIonly$CellType)

markersListP.str<-split(markersP.str,markersP.str$CellType)
markersListP.str.BrainOnly<-split(markersP.str.BrainOnly,markersP.str.BrainOnly$CellType)


markersC<-read.delim("CellMarkerMouse_cell_markers.txt")
markersC<-markersC[markersC$cancerType=="Normal",]
markersC$tissueType[markersC$tissueType=="Gastrointestinal tract"]<-"GI tract"
markersC$CellType<-paste0(markersC$cellName," (",markersC$tissueType,")")
markersC$db<-"CellMarker"
unique(markersC$CellType)
exclude<-grepl("Schwalie",markersC$CellType)
markersC<-markersC[!exclude,]


C.tissues.dc.all<-c("GI tract","Intestinal crypt",
                    "Colon epithelium",
                    "Ileum","Blood","Peyer patch","Neural tube",
                    "Intestine","Lymphoid Tissue",
                    "Colon","Undefined",
                    "Bone Marrow","Adipose tissue",
                    "Peritoneal cavity",
                    "Spinal cord","Lymph node",
                    "Blood vessel","Artery",
                    "Epithelium","Muscle" ,
                    "Small intestine" ,"Brain",
                    "Mesenteric lymph node")
C.tissues.dc.gi<-c("GI tract","Intestinal crypt",
                   "Colon epithelium",
                   "Ileum","Peyer patch",
                   "Intestine","Colon","Small intestine")

C.tissues.str.all<-c("Brain","Blood","Blood vessel","Artery",
                     "Undefined","Neural tube")
C.tissues.str.brain<-c("Brain")

markersC.dc<-markersC[markersC$tissueType%in%C.tissues.dc.all,]


markersC.dc.GIonly<-markersC[markersC$tissueType%in%C.tissues.dc.gi,]

markersC.str<-markersC[markersC$tissueType%in%C.tissues.str.all,]

markersC.str.BrainOnly<-markersC[markersC$tissueType%in%C.tissues.str.brain,]

markersListC.dc<-split(markersC.dc,markersC.dc$CellType)
markersListC.dc.symbs<-lapply(markersListC.dc,function(x){
  symbs<-unlist(strsplit(x$geneSymbol,", "))
  return(symbs)
  
})


markersListC.dc.GIonly<-split(markersC.dc.GIonly,markersC.dc.GIonly$CellType)
markersListC.dc.symbs.GIonly<-lapply(markersListC.dc.GIonly,function(x){
  symbs<-unlist(strsplit(x$geneSymbol,", "))
  return(symbs)
  
})

markersListC.str<-split(markersC.str,markersC.str$CellType)
markersListC.str.symbs<-lapply(markersListC.str,function(x){
  symbs<-unlist(strsplit(x$geneSymbol,", "))
  return(symbs)
  
})

markersListC.str.BrainOnly<-split(markersC.str.BrainOnly,markersC.str.BrainOnly$CellType)
markersListC.str.symbs.BrainOnly<-lapply(markersListC.str.BrainOnly,function(x){
  symbs<-unlist(strsplit(x$geneSymbol,", "))
  return(symbs)
  
})

# dc-------
mm=0.2
tiss="dc"
markerSet<-list(GI.only=list(list.P=markersListP.dc.GIonly,
                             list.C=markersListC.dc.GIonly,
                             list.C.symbs=markersListC.dc.symbs.GIonly),
                All.relevant=list(list.P=markersListP.dc,
                                  list.C=markersListC.dc,
                                  list.C.symbs=markersListC.dc.symbs)
)

gi<-readRDS(file.path(tiss,"geneINFOall.rds"))
gi<-gi[!is.na(gi$geneSymbol),]
allgenes<-gi$geneSymbol
n.all<-length(allgenes)
gi<-gi[gi$MM>=mm,]
mods<-as.character(readRDS(file.path(tiss,"modtab.rds"))$moduleColor)


overlapTableP<-NULL
overlapTableC<-NULL

for(markerset in names(markerSet)){
  list.P<-markerSet[[markerset]][["list.P"]]
  list.C<-markerSet[[markerset]][["list.C"]]
  list.C.symbs<-markerSet[[markerset]][["list.C.symbs"]]
  overlapTableP0<-NULL
  overlapTableC0<-NULL
  for(i in 1:length(list.P)){
    for(j in 1:length(mods)){
      cellType<-names(list.P)[i]
      cellgenes<-list.P[[cellType]]$mgi_symbol
      cellgenes.inallgenes<-intersect(cellgenes,allgenes)
      
      module<-mods[j]
      modgenes<-gi$geneSymbol[gi$moduleColor==module]
      
      n.module<-length(modgenes)
      n.geneset<-length(cellgenes.inallgenes)
      overlap<-intersect(cellgenes.inallgenes,modgenes)
      n.overlap<-length(overlap)
      
      pval<-phyper(n.overlap-1, n.module, n.all-n.module, 
                   n.geneset,lower.tail= FALSE)
      
      ovlerlapSymbs<-gi$geneSymbol[gi$geneSymbol%in%overlap]
      
      overlapTableP0<-rbind(overlapTableP0,
                            c(list.P[[cellType]]$organ[1],
                              list.P[[cellType]]$cell.type[1],
                              cellType,module,n.overlap,pval,n.geneset,n.module,
                              paste(ovlerlapSymbs,collapse = "/"),
                              list.P[[cellType]]$db[1],
                              gsub("\\."," ",markerset)))
    }
  }
  
  for(i in 1:length(list.C)){
    for(j in 1:length(mods)){
      cellType<-names(list.C)[i]
      cellgenes<-list.C.symbs[[cellType]]
      cellgenes.inallgenes<-intersect(cellgenes,allgenes)
      
      module<-mods[j]
      modgenes<-gi$geneSymbol[gi$moduleColor==module]
      
      n.module<-length(modgenes)
      n.geneset<-length(cellgenes.inallgenes)
      overlap<-intersect(cellgenes.inallgenes,modgenes)
      n.overlap<-length(overlap)
      pval<-phyper(n.overlap-1, n.module, n.all-n.module, 
                   n.geneset,lower.tail= FALSE)
      
      ovlerlapSymbs<-gi$geneSymbol[gi$geneSymbol%in%overlap]
      overlapTableC0<-rbind(overlapTableC0,
                            c(list.C[[cellType]]$tissueType[1],
                              list.C[[cellType]]$cellName[1],
                              cellType,module,n.overlap,pval,n.geneset,n.module,
                              paste(ovlerlapSymbs,collapse = "/"),
                              list.C[[cellType]]$db[1],
                              gsub("\\."," ",markerset)))
    }
  }
  
  overlapTableP0<-as.data.frame(overlapTableP0)
  names(overlapTableP0)<-c("Tissue","CellType","CellTypeFull","Module","Overlap","P","SetSize","ModSize","CommonGenes","Database","MarkerSet")
  
  overlapTableP0$P<-as.numeric(overlapTableP0$P)
  overlapTableP0<-overlapTableP0[overlapTableP0$Overlap>0,]
  for(mod in mods){
    overlapTableP0[overlapTableP0$Module==mod,"p.adj"]<-p.adjust(overlapTableP0[overlapTableP0$Module==mod,"P"],method = "BH")
  }
  
  
  overlapTableC0<-as.data.frame(overlapTableC0)
  names(overlapTableC0)<-c("Tissue","CellType","CellTypeFull","Module","Overlap","P","SetSize","ModSize","CommonGenes","Database","MarkerSet")
  
  overlapTableC0$P<-as.numeric(overlapTableC0$P)
  overlapTableC0<-overlapTableC0[overlapTableC0$Overlap>0,]
  for(mod in mods){
    overlapTableC0[overlapTableC0$Module==mod,"p.adj"]<-p.adjust(overlapTableC0[overlapTableC0$Module==mod,"P"],method = "BH")
  }
  
  if(markerset=="All.relevant"){
    overlapTableP<-rbind(overlapTableP,
                         overlapTableP0[!overlapTableP0$Tissue%in%P.tissues.dc.gi,])
    overlapTableC<-rbind(overlapTableC,
                         overlapTableC0[!overlapTableC0$Tissue%in%C.tissues.dc.gi,])
    
  }else{
    overlapTableP<-rbind(overlapTableP,overlapTableP0)
    overlapTableC<-rbind(overlapTableC,overlapTableC0)
  }
  
}

overlapTable.dc<-rbind(overlapTableP,overlapTableC)
overlapTable.dc<-overlapTable.dc[order(overlapTable.dc$Module,overlapTable.dc$P),]
write.csv(overlapTable.dc,row.names = F,"DcCellTypeHyperG.csv")



# str --------
tiss="str"

gi<-readRDS(file.path(tiss,"geneINFOall.rds"))
gi<-gi[!is.na(gi$geneSymbol),]
allgenes<-gi$geneSymbol
n.all<-length(allgenes)
gi<-gi[gi$MM>=mm,]
mods<-as.character(readRDS(file.path(tiss,"modtab.rds"))$moduleColor)



markerSet<-list(Brain.only=list(list.P=markersListP.str.BrainOnly,
                                list.C=markersListC.str.BrainOnly,
                                list.C.symbs=markersListC.str.symbs.BrainOnly),
                All.relevant=list(list.P=markersListP.str,
                                  list.C=markersListC.str,
                                  list.C.symbs=markersListC.str.symbs)
)


overlapTableP<-NULL
overlapTableC<-NULL

for(markerset in names(markerSet)){
  list.P<-markerSet[[markerset]][["list.P"]]
  list.C<-markerSet[[markerset]][["list.C"]]
  list.C.symbs<-markerSet[[markerset]][["list.C.symbs"]]
  overlapTableP0<-NULL
  overlapTableC0<-NULL
  for(i in 1:length(list.P)){
    for(j in 1:length(mods)){
      cellType<-names(list.P)[i]
      cellgenes<-list.P[[cellType]]$mgi_symbol
      cellgenes.inallgenes<-intersect(cellgenes,allgenes)
      
      module<-mods[j]
      modgenes<-gi$geneSymbol[gi$moduleColor==module]
      
      n.module<-length(modgenes)
      n.geneset<-length(cellgenes.inallgenes)
      overlap<-intersect(cellgenes.inallgenes,modgenes)
      n.overlap<-length(overlap)
      
      pval<-phyper(n.overlap-1, n.module, n.all-n.module, 
                   n.geneset,lower.tail= FALSE)
      
      ovlerlapSymbs<-gi$geneSymbol[gi$geneSymbol%in%overlap]
      
      overlapTableP0<-rbind(overlapTableP0,
                            c(list.P[[cellType]]$organ[1],
                              list.P[[cellType]]$cell.type[1],
                              cellType,module,n.overlap,pval,n.geneset,n.module,
                              paste(ovlerlapSymbs,collapse = "/"),
                              list.P[[cellType]]$db[1],
                              gsub("\\."," ",markerset)))
    }
  }
  
  for(i in 1:length(list.C)){
    for(j in 1:length(mods)){
      cellType<-names(list.C)[i]
      cellgenes<-list.C.symbs[[cellType]]
      cellgenes.inallgenes<-intersect(cellgenes,allgenes)
      
      module<-mods[j]
      modgenes<-gi$geneSymbol[gi$moduleColor==module]
      
      n.module<-length(modgenes)
      n.geneset<-length(cellgenes.inallgenes)
      overlap<-intersect(cellgenes.inallgenes,modgenes)
      n.overlap<-length(overlap)
      pval<-phyper(n.overlap-1, n.module, n.all-n.module, 
                   n.geneset,lower.tail= FALSE)
      
      ovlerlapSymbs<-gi$geneSymbol[gi$geneSymbol%in%overlap]
      overlapTableC0<-rbind(overlapTableC0,
                            c(list.C[[cellType]]$tissueType[1],
                              list.C[[cellType]]$cellName[1],
                              cellType,module,n.overlap,pval,n.geneset,n.module,
                              paste(ovlerlapSymbs,collapse = "/"),
                              list.C[[cellType]]$db[1],
                              gsub("\\."," ",markerset)))
    }
  }
  
  overlapTableP0<-as.data.frame(overlapTableP0)
  names(overlapTableP0)<-c("Tissue","CellType","CellTypeFull","Module","Overlap","P","SetSize","ModSize","CommonGenes","Database","MarkerSet")
  
  overlapTableP0$P<-as.numeric(overlapTableP0$P)
  overlapTableP0<-overlapTableP0[overlapTableP0$Overlap>0,]
  for(mod in mods){
    overlapTableP0[overlapTableP0$Module==mod,"p.adj"]<-p.adjust(overlapTableP0[overlapTableP0$Module==mod,"P"],method = "BH")
  }
  
  
  overlapTableC0<-as.data.frame(overlapTableC0)
  names(overlapTableC0)<-c("Tissue","CellType","CellTypeFull","Module","Overlap","P","SetSize","ModSize","CommonGenes","Database","MarkerSet")
  
  overlapTableC0$P<-as.numeric(overlapTableC0$P)
  overlapTableC0<-overlapTableC0[overlapTableC0$Overlap>0,]
  for(mod in mods){
    overlapTableC0[overlapTableC0$Module==mod,"p.adj"]<-p.adjust(overlapTableC0[overlapTableC0$Module==mod,"P"],method = "BH")
  }
  
  if(markerset=="All.relevant"){
    overlapTableP<-rbind(overlapTableP,
                         overlapTableP0[!overlapTableP0$Tissue%in%P.tissues.str.brain,])
    overlapTableC<-rbind(overlapTableC,
                         overlapTableC0[!overlapTableC0$Tissue%in%C.tissues.str.brain,])
    
  }else{
    overlapTableP<-rbind(overlapTableP,overlapTableP0)
    overlapTableC<-rbind(overlapTableC,overlapTableC0)
  }
  
}

overlapTable.str<-rbind(overlapTableP,overlapTableC)
overlapTable.str<-overlapTable.str[order(overlapTable.str$Module,overlapTable.str$P),]
write.csv(overlapTable.str,row.names = F,"strCellTypeHyperG.csv")
