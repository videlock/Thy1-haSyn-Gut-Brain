#  Resources requested:
#$ -l h_data=2G,h_rt=00:30:00
#$ -V
#  Email address to notify
#$ -M videlock@mail
#$ -m a


rm(list=ls());
suppressMessages(library(GOstats,quietly = T))
suppressMessages(library(org.Mm.eg.db,quietly = T))

# options(echo=TRUE)
options(stringsAsFactors = FALSE);

onts<-c("BP","CC","MF")
MMthreshs<-c("6","7","8","all")




mods<-as.character(readRDS("modtab.rds")$moduleColor)


param.df<-data.frame()

for(mod in mods){
  for(ont in onts){
    for(MMthresh in MMthreshs){
      param.df<-rbind(param.df,
                      c(mod,ont,MMthresh))
    }
  }
}


names(param.df)<-c("mod","ont","MMthresh")

taskID = as.integer(Sys.getenv("SGE_TASK_ID"))


ont=param.df$ont[taskID]
MMthresh=param.df$MMthresh[taskID]
mod=param.df$mod[taskID]

print(ont)
print(MMthresh)
print(mod)

load("processed_data/inputData.rda")


gi<-readRDS("geneINFO_MM6.rds")

genedata<-an
allLLIDs <- genedata$NCBI.gene.ID
allLLIDs<-allLLIDs[!allLLIDs%in%NA]
allLLIDs<-as.character(allLLIDs)


gi<-gi[gi$moduleColor==mod,]

go.outMain<-"GO"
dir.create(go.outMain,showWarnings = F)

go.out.ont<-file.path(go.outMain,ont)
dir.create(go.out.ont,showWarnings = F)

go.out.ont.mm<-file.path(go.out.ont,paste0("MM",MMthresh))
dir.create(go.out.ont.mm,showWarnings = F)


if(MMthresh=="all"){
  modgenes<-gi$EntrezID
}else{
  gi2<-gi[gi$MM>=as.numeric(MMthresh)/10,]
  modgenes<-gi2$EntrezID
}

modgenes<-modgenes[!modgenes%in%NA]
modgenes<-as.character(modgenes)

if(length(modgenes)>20){
  print("running GO")
  
  GOparams <- new("GOHyperGParams", 
                  geneIds = modgenes, 
                  universeGeneIds = allLLIDs, 
                  annotation = "org.Mm.eg.db", 
                  ontology = ont, pvalueCutoff = 0.05,
                  conditional = T, testDirection = "over")
  hgOver <- hyperGTest(GOparams)
  gotable <- summary(hgOver)
  
  saveRDS(gotable,file = file.path(go.out.ont.mm,paste0(mod,".rds")))
  
}
