library(tidyverse)
gi.dc <-
  readRDS("~/Dropbox/RNAseq/ASO_QuantSeq_Final/WGCNA/dc/geneINFOall.rds")


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

# names(gi.gries)[22]<-"fc.li"

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

# s1m<-m1m&gi.gries$padj1m<0.1
# s3m<-m3m&gi.gries$padj3m<0.1
# 
# sig<-ifelse(s1m&s3m,
#             "1m,3m",
#             ifelse(s1m,
#                    "1m",
#                    ifelse(s3m,
#                           "3m",
#                           NA)))
# 
# 
# gi.gries$Sig<-sig
# 
# 
# s1mp<-m1m&gi.gries$P.Value_1m<0.05
# s3mp<-m3m&gi.gries$P.Value_3m<0.05
# 
# sigPmatched<-ifelse(s1mp&s3mp,
#              "1m,3m",
#              ifelse(s1mp,
#                     "1m",
#                     ifelse(s3mp,
#                            "3m",
#                            NA)))
# 
# 
# gi.gries$SigNomPmatched<-sigPmatched
# 
# gi.gries$padj1mAll<-p.adjust(gi.gries$P.Value_1m, method = "BH")
# 
# gi.gries$padj3mAll<-p.adjust(gi.gries$P.Value_3m, method = "BH")
# 
# gi.gries$padj3m[gi.gries$MatchDir%in%c("3m","1m,3m")]<-p.adjust(gi.gries$P.Value_3m[gi.gries$MatchDir%in%c("3m","1m,3m")],method = "BH")

# 
# s1mpall<-gi.gries$padj1mAll<0.1
# s3mpall<-gi.gries$padj3mAll<0.1
# 
# sigPAll<-ifelse(s1mpall&s3mpall,
#                     "1m,3m",
#                     ifelse(s1mpall,
#                            "1m",
#                            ifelse(s3mpall,
#                                   "3m",
#                                   NA)))
# 
# 
# gi.gries$SigAny<-sigPAll
# 
# s1mpallp<-gi.gries$P.Value_1m<0.05
# s3mpallp<-gi.gries$P.Value_3m<0.05
# 
# sigPAllNom<-ifelse(s1mpallp&s3mpallp,
#                 "1m,3m",
#                 ifelse(s1mpallp,
#                        "1m",
#                        ifelse(s3mpallp,
#                               "3m",
#                               NA)))
# 
# 
# gi.gries$SigAnyNomP<-sigPAllNom
# 
# 
# for(j in 1:ncol(gi.gries)){
#   if(is.numeric(gi.gries[,j])){
#     gi.gries[,j]<-round(gi.gries[,j],digits = 2)
#   }
# }

# names(gi.gries)

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


saveRDS(gi.gries,file="~/Dropbox/Papers/ASO_BrainGutGeneExpression/PaperRmd/Tables/griesmatching.rds")
