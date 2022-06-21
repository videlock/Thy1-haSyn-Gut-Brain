
options(suppressPackageStartupMessages=T)
library(Biobase)
library(biomaRt)

load("../WGCNA/WGCNA_Apr2021/cons/inputData.rda")
dc.expr<-t(scale(multiExpr$dc$data))
str.expr<-t(scale(multiExpr$str$data))

metadata.dc<-metadataList$dc
rownames(metadata.dc)<-paste0("dc_",rownames(metadata.dc))
colnames(dc.expr)<-rownames(metadata.dc)

metadata.str<-metadataList$str
rownames(metadata.str)<-paste0("str_",rownames(metadata.str))
colnames(str.expr)<-rownames(metadata.str)

datexpr<-cbind(dc.expr,str.expr)
metadata<-rbind(metadata.dc,metadata.str)
gt<-ifelse(metadata$GT=="Hem","A","W")
tissue<-ifelse(metadata$tissue=="dc","C","S")
grp<-paste(gt,metadata$Time,tissue,sep = "")
# rownames(metadata)
# colnames(datexpr)


sink(file = file.path("grp.cls"))
cat(paste(length(grp), length(levels(as.factor(grp))), 1))
cat("\n")
cat("#",paste(unique(grp), collapse = " "))
cat("\n")
cat(grp)
sink()

musGenes <- rownames(datexpr)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
annot = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = musGenes , mart = mouse, attributesL = c("hgnc_symbol","description"), martL = human, uniqueRows=T)
annot<-annot[!annot$HGNC.symbol=="",]
annot$Gene.description<-gsub("[[:space:]]\\[.+","",annot$Gene.description)

# mouse ids linked to multiple human
ambigIDs<-unique(annot$Gene.stable.ID[duplicated(annot$Gene.stable.ID)]) #121
length(unique(annot$HGNC.symbol)) 
ambigSymbs<-unique(annot$HGNC.symbol[duplicated(annot$HGNC.symbol)]) #21
isambig<-annot$Gene.stable.ID%in%ambigIDs|annot$HGNC.symbol%in%ambigSymbs
ambig<-annot[isambig,]
ambig$isDupSymb<-ifelse(ambig$HGNC.symbol%in%ambigSymbs,"dup","")
write.csv(ambig[order(ambig$Gene.stable.ID),],row.names = F,"amigIDsAndSymbs.csv")
ambig<-read.csv("amigIDsAndSymbs.csv")

annot1<-annot[!isambig,]
annot2<-ambig[ambig$keep==1&!is.na(ambig$keep),names(annot1)]
annot3<-ambig[is.na(ambig$keep),names(annot1)]
annotNew<-rbind(annot1,annot2,annot3)
inannotNew<-rownames(datexpr)%in%annotNew$Gene.stable.ID
datexpr2<-datexpr[inannotNew,]
dim(datexpr2)
rownames(annotNew)<-annotNew$Gene.stable.ID
annotNew2<-annotNew[rownames(datexpr2),]
library(WGCNA)
Col.datexpr2<-collapseRows(datexpr2,
                           rowGroup = annotNew$HGNC.symbol,
                           rowID = annotNew$Gene.stable.ID)

annot.final<-annotNew2[Col.datexpr2$selectedRow,]
datexpr.final<-datexpr2[Col.datexpr2$selectedRow,]
datexpr.symbs<-datexpr.final
rownames(datexpr.symbs)<-annot.final$HGNC.symbol
mman<-an[rownames(annot.final),]
mman$HSsymb<-annot.final$HGNC.symbol
write.csv(mman,"annotations.csv",row.names = F)

write.table(annot.final, file = file.path("fdat.chip"), row.names = F,
            sep = "\t", col.names = c("Probe Set ID","Gene Symbol","Gene Title"),
            quote = F)       



frontcols<-data.frame(NAME=rownames(datexpr.symbs),DESCRIPTION=NA)
expdat<-cbind(frontcols,datexpr.symbs)
write.table(expdat, sep = "\t", quote = F, 
            file = file.path("exprs.txt"), row.names = F)



pathBase="~/project-videlock/QuantSeqApr2019/FullWorkflow/E_GSEA"

genesets=c(hallmark="h.all.v7.0.symbols.gmt",
           GOBP="c5.bp.v7.0.symbols.gmt",
           KEGG="c2.cp.kegg.v7.0.symbols.gmt")

genesetPath="~/project-videlock/GSEA/GeneSets"


# cmdFileDir=paste0(pathBase,"runFilesGTtime/")
cmdFileDir="runFiles"

dir.create(cmdFileDir,showWarnings = FALSE)

comps=c("A1C_versus_W1C",
        "A3C_versus_W3C",
        "A1S_versus_W1S",
        "A3S_versus_W3S",
        "A1C_versus_A1S",
        "A3C_versus_A3S"
)

gseaCmds<-character()
Labels<-character()

inPathMain=pathBase

for(i in 1:length(comps)){
  for(j in 1:length(genesets)){
    comp=comps[i]
    gmxPath=file.path(genesetPath,genesets[j])
    rptLabel=paste(comp,names(genesets)[j],sep = "_")
    gseaCmd<-paste("./gsea-cli.sh GSEA",
                   "-res",
                   file.path(inPathMain,"exprs.txt"),
                   "-cls",
                   file.path(inPathMain,paste0("grp.cls#",comp)),
                   paste0(" -gmx ",gmxPath),
                   "-collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted",
                   paste0("-rpt_label ", rptLabel),
                   "-metric Signal2Noise -sort real -order descending",
                   # "-chip ",
                   # file.path(inPathMain,"fdat.chip"),
                   "-create_gcts true -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists tru -set_max 500 -set_min 15 -zip_report true",
                   paste0("-out ",inPathMain)
    )
                   
                   gseaCmds<-c(gseaCmds,gseaCmd)
                   Labels<-c(Labels,rptLabel)
                   
  }
}



for(i in 1:length(gseaCmds)){
  sink(file = file.path(cmdFileDir,paste0(Labels[i],"_",i,".sh")))
  cat("#! bin/bash")
  cat("\n\n")
  cat(gseaCmds[i])
  sink()
}


write.table(Labels,"grp_cmdref.txt",quote = F)

sink(file = "runGSEA.sh")
cat("#!/bin/bash")
cat("\n\n")
cat("#$ -cwd\n")
cat("#$ -o  logs/runGsea$JOB_ID_$TASK_ID.txt\n")
cat("#$ -j y\n")
cat("#$ -t 1:",length(gseaCmds),"\n",sep = "")
cat("#$ -l h_data=16G,h_rt=01:00:00\n")
cat("#$ -V\n")
cat("#$ -M videlock@mail\n")
cat("#$ -m a\n")
cat("\n\n")
cat("\n\n") 
cat(". /u/local/Modules/default/init/modules.sh")
cat("\n\n")
cat("module load java/13.0.1")
cat("\n\n")
cat("cd /u/project/videlock/videlock/GSEA/GSEA_Linux_4.0.3 && source  /u/project/videlock/videlock/QuantSeqApr2019/FullWorkflow/E_GSEA/runFiles/*_${SGE_TASK_ID}.sh")
sink()




comps2=c(
        "W1C_versus_W1S",
        "W3C_versus_W3S",
        "A3C_versus_A1C",
        "W3C_versus_W1C",
        "A3S_versus_A1S",
        "W3S_versus_W1S"
)

gseaCmds<-character()
Labels<-character()

inPathMain=pathBase

for(i in 1:length(comps)){
  for(j in 1:length(genesets)){
    comp=comps[i]
    gmxPath=file.path(genesetPath,genesets[j])
    rptLabel=paste(comp,names(genesets)[j],sep = "_")
    gseaCmd<-paste("./gsea-cli.sh GSEA",
                   "-res",
                   file.path(inPathMain,"exprs.txt"),
                   "-cls",
                   file.path(inPathMain,paste0("grp.cls#",comp)),
                   paste0(" -gmx ",gmxPath),
                   "-collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted",
                   paste0("-rpt_label ", rptLabel),
                   "-metric Signal2Noise -sort real -order descending",
                   # "-chip ",
                   # file.path(inPathMain,"fdat.chip"),
                   "-create_gcts true -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists tru -set_max 500 -set_min 15 -zip_report true",
                   paste0("-out ",inPathMain)
    )
    
    gseaCmds<-c(gseaCmds,gseaCmd)
    Labels<-c(Labels,rptLabel)
    
  }
}

cmdFileDir="runFiles2"

dir.create(cmdFileDir,showWarnings = FALSE)

for(i in 1:length(gseaCmds)){
  sink(file = file.path(cmdFileDir,paste0(Labels[i],"_",i,".sh")))
  cat("#! bin/bash")
  cat("\n\n")
  cat(gseaCmds[i])
  sink()
}


write.table(Labels,"grp_cmdref2.txt",quote = F)

sink(file = "runGSEA2.sh")
cat("#!/bin/bash")
cat("\n\n")
cat("#$ -cwd\n")
cat("#$ -o  logs2/runGsea$JOB_ID_$TASK_ID.txt\n")
cat("#$ -j y\n")
cat("#$ -t 1:",length(gseaCmds),"\n",sep = "")
cat("#$ -l h_data=16G,h_rt=01:00:00\n")
cat("#$ -V\n")
cat("#$ -M videlock@mail\n")
cat("#$ -m a\n")
cat("\n\n")
cat("\n\n") 
cat(". /u/local/Modules/default/init/modules.sh")
cat("\n\n")
cat("module load java/13.0.1")
cat("\n\n")
cat("cd /u/project/videlock/videlock/GSEA/GSEA_Linux_4.0.3 && source  /u/project/videlock/videlock/QuantSeqApr2019/FullWorkflow/E_GSEA/runFiles2/*_${SGE_TASK_ID}.sh")
sink()


# all--------
dirmain="all"
dir.create(dirmain)

sink(file = file.path(dirmain,"grp.cls"))
cat(paste(length(grp), length(levels(as.factor(grp))), 1))
cat("\n")
cat("#",paste(unique(grp), collapse = " "))
cat("\n")
cat(grp)
sink()

write.table(expdat, sep = "\t", quote = F, 
            file = file.path(dirmain,"exprs.txt"), row.names = F)


pathBase="~/project-videlock/QuantSeqApr2019/FullWorkflow/E_GSEA/all"
genesetsAll<-c("h.all.v7.0.symbols.gmt","c5.bp.v7.0.symbols.gmt","c2.cp.kegg.v7.0.symbols.gmt")
genesetPath="~/project-videlock/GSEA/GeneSets"
cmdFileDir=file.path(dirmain,"runFiles")
dir.create(cmdFileDir,showWarnings = FALSE)


comps=c("A1C_versus_W1C",
        "A3C_versus_W3C",
        "A1S_versus_W1S",
        "A3S_versus_W3S",
        "A1C_versus_A1S",
        "A3C_versus_A3S",
        "W1C_versus_W1S",
        "W3C_versus_W3S",
        "A3C_versus_A1C",
        "W3C_versus_W1C",
        "A3S_versus_A1S",
        "W3S_versus_W1S")

gseaCmds<-character()
Labels<-character()

inPathMain=pathBase
genesetPath="/u/project/videlock/videlock/GSEA/GeneSets"
gmxPath=paste(file.path(genesetPath,genesetsAll),collapse = ",")


for(i in 1:length(comps)){
  comp=comps[i]
  rptLabel=comp
  gseaCmd<-paste("./gsea-cli.sh GSEA",
                 "-res",
                 file.path(inPathMain,"exprs.txt"),
                 "-cls",
                 file.path(inPathMain,paste0("grp.cls#",comp)),
                 paste0(" -gmx ",gmxPath),
                 "-collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted",
                 paste0("-rpt_label ", rptLabel),
                 "-metric Signal2Noise -sort real -order descending",
                 # "-chip ",
                 # file.path(inPathMain,"fdat.chip"),
                 "-create_gcts true -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists tru -set_max 500 -set_min 15 -zip_report true",
                 paste0("-out ",inPathMain)
  )
  
  gseaCmds<-c(gseaCmds,gseaCmd)
  Labels<-c(Labels,rptLabel)
}



for(i in 1:length(gseaCmds)){
  sink(file = file.path(cmdFileDir,paste0(Labels[i],"_",i,".sh")))
  cat("#! bin/bash")
  cat("\n\n")
  cat(gseaCmds[i])
  sink()
}


write.table(Labels,file.path(dirmain,"grp_cmdref.txt"),quote = F)

sink(file = file.path(dirmain,"runGSEA.sh"))
cat("#!/bin/bash")
cat("\n\n")
cat("#$ -cwd\n")
cat("#$ -o  logs/runGsea$JOB_ID_$TASK_ID.txt\n")
cat("#$ -j y\n")
cat("#$ -t 1:",length(gseaCmds),"\n",sep = "")
cat("#$ -l h_data=16G,h_rt=01:00:00\n")
cat("#$ -V\n")
cat("#$ -M videlock@mail\n")
cat("#$ -m a\n")
cat("\n\n")
cat("\n\n") 
cat(". /u/local/Modules/default/init/modules.sh")
cat("\n\n")
cat("module load java/13.0.1")
cat("\n\n")
cat("cd /u/project/videlock/videlock/GSEA/GSEA_Linux_4.0.3 && source  /u/project/videlock/videlock/QuantSeqApr2019/FullWorkflow/E_GSEA/all/runFiles/*_${SGE_TASK_ID}.sh")
sink()


