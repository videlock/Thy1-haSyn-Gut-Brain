library(tidyverse)
gi.dc <-
  readRDS("../WGCNA/dc/geneINFOall.rds")


gries<-read.csv("griesprot.csv")


keep=gries$p.li<0.05
gries.li<-gries[keep,]

protnames<-strsplit(gries.li$Protein.Name,";")

uniprot_mapping <- function(ids) {
  uri <- 'http://www.uniprot.org/uniprot/?query='
  idStr <- paste(ids, collapse="+or+")
  format <- '&format=tab'
  fullUri <- paste0(uri,idStr,format)
  dat <- read.delim(fullUri)
  dat
}

protannot<-list()
for(i in 1:length(protnames)){
  pids<-protnames[[i]]
  updat<-uniprot_mapping(pids)
  protannot[[i]]<-updat
}

for(i in 1:length(protannot)){
  df<-protannot[[i]]
  df$rownum<-i
  protannot[[i]]<-df
}

protannotTab<-as.data.frame(do.call(rbind,protannot))
write.csv(protannotTab,"protannotTab.csv")

protannotTab<-read.csv("protannotTab.csv")

gries.li$geneSymbol<-protannotTab$Gene.names

rownames(gi.dc) <- NULL

gi.gries<-inner_join(gi.dc,gries.li,"geneSymbol")


m1m<-(gi.gries$logFC_1m<0 &
            gi.gries$fc.li<0)|
  (gi.gries$logFC_1m>0 &
     gi.gries$fc.li>0)

m3m<-(gi.gries$logFC_3m<0 &
            gi.gries$fc.li<0)|
  (gi.gries$logFC_3m>0 &
     gi.gries$fc.li>0)


matchdir<-ifelse(m1m&m3m,
                "Yes (1m,3m)",
                ifelse(m1m,
                       "Yes (1m)",
                       ifelse(m3m,
                              "Yes (3m)",
                              "No")))

gi.gries$MatchDir<-matchdir

gi.gries$padj1m[!gi.gries$MatchDir%in%c("No")]<-p.adjust(gi.gries$P.Value_1m[!gi.gries$MatchDir%in%c("No")],method = "BH")

gi.gries$padj3m[!gi.gries$MatchDir%in%c("No")]<-p.adjust(gi.gries$P.Value_3m[!gi.gries$MatchDir%in%c("No")],method = "BH")

write.csv(gi.gries[,c("geneSymbol", "Name",
                      "MatchDir",

                      "fc.li","p.li",
                      "logFC_1m", 
                      "P.Value_1m",
                      "padj1m" ,
                      "logFC_3m" ,
 
                      "P.Value_3m",
                      "padj3m",
                      "moduleColor","MM",
                      "Protein.Name",
                      "ENSG" 

                      )],
          "griesmatch.csv",row.names = F,na = "")

