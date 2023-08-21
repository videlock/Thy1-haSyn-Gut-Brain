

options(echo=TRUE)
options(stringsAsFactors = FALSE);
library(WGCNA);
allowWGCNAThreads() 
library(flashClust);
library(gplots)
library(RColorBrewer)

setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))

# striatum-----

outdir="CSstr"
dir.create(outdir)

load("Striatum/processed_data/inputData.rda")

#params------
maxPout=1
sp=12
ch<-0.995
mms<-60
ds<-1
dthresh<-0.25


adjacency = adjacency(datExpr, power = sp, type = "signed",
                      corFnc="bicor",corOptions = list(maxPOutliers =maxPout/100))
TOM = TOMsimilarity(adjacency);

gi<-readRDS("str/geneINFO_MM6.rds")
gi<-gi[!is.na(gi$EntrezID),]
gi$color<-col2hex(gi$moduleColor)

genes = names(datExpr)
selected = genes%in%gi$ENSG
selectedGenes = genes[selected]
gi<-gi[match(selectedGenes,gi$ENSG),]
selectedTOM = TOM[selected, selected];
dimnames(selectedTOM) = list(selectedGenes, selectedGenes)

selectedAdj=adjacency[selected, selected];
dimnames(selectedAdj) = list(selectedGenes, selectedGenes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(selectedAdj,
                               edgeFile = file.path(outdir,"Edges.txt"),
                               nodeFile = file.path(outdir,"Nodes.txt"),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = selectedGenes,
                               altNodeNames = gi$geneSymbol,
                               nodeAttr = gi);



write.table(as.data.frame(an$Gene.name),file.path(outdir,"background.txt"),row.names = F,quote = F,col.names = F)


# colon --------

outdir="CSdc"
dir.create(outdir)

load("dc/processed_data/inputData.rda")

#params------
maxPout=5
sp=8
ch<-0.995
mms<-60
ds<-1
dthresh<-0.25

adjacency = adjacency(datExpr, power = sp, type = "signed",
                      corFnc="bicor",corOptions = list(maxPOutliers =maxPout/100))

gi<-readRDS("dc/geneINFO_MM6.rds")
gi<-gi[!is.na(gi$EntrezID),]
gi$color<-col2hex(gi$moduleColor)

genes = names(datExpr)
selected = genes%in%gi$ENSG
selectedGenes = genes[selected]
gi<-gi[match(selectedGenes,gi$ENSG),]

selectedAdj=adjacency[selected, selected];
dimnames(selectedAdj) = list(selectedGenes, selectedGenes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(selectedAdj,
                               edgeFile = file.path(outdir,"Edges-dc-mm6.txt"),
                               nodeFile = file.path(outdir,"Nodes-dc-mm6.txt"),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = selectedGenes,
                               altNodeNames = gi$geneSymbol,
                               nodeAttr = gi);


modColortable<-data.frame(sharedName=paste(cyt$edgeData$fromNode,
                                           " (interacts with) ",
                                           cyt$edgeData$toNode),
                          fromColor=cyt$edgeData$fromAltName,
                          toColor=cyt$edgeData$toAltName)
write.csv(modColortable,file.path(outdir,"modcolortable.csv"),row.names = F)



write.table(as.data.frame(an$Gene.name),file.path(outdir,"background.txt"),row.names = F,quote = F,col.names = F)





# consensus-------

outdir="CScons"
dir.create(outdir)

load("cons/inputData.rda")

write.table(as.data.frame(an$Gene.name),file.path(outdir,"background.txt"),row.names = F,quote = F,col.names = F)

adjacency.dc = adjacency(multiExpr$dc$data, power = 8, type = "signed",
                         corFnc="bicor",corOptions = list(maxPOutliers =0.05))

adjacency.str = adjacency(multiExpr$str$data, power = 12, type = "signed",
                          corFnc="bicor",corOptions = list(maxPOutliers =0.01))


gi<-readRDS("cons/geneINFOallgenes.rds")

gi$color<-col2hex(gi$moduleColor)

gi.dc<-gi[gi$MM.dc>0.6,]
gi.str<-gi[gi$MM.str>0.6,]

genes = names(multiExpr$dc$data)
selected.dc = genes%in%gi.dc$ENSG
selectedGenes.dc = genes[selected.dc]
gi.dc<-gi.dc[match(selectedGenes.dc,gi.dc$ENSG),]

selectedAdj.dc=adjacency.dc[selected.dc, selected.dc];
dimnames(selectedAdj.dc) = list(selectedGenes.dc, selectedGenes.dc)


cyt = exportNetworkToCytoscape(selectedAdj.dc,
                               edgeFile = file.path(outdir,"Edges_ConsAdjDc.txt"),
                               nodeFile = file.path(outdir,"Nodes_ConsAdjDc.txt"),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = selectedGenes.dc,
                               altNodeNames = gi.dc$geneSymbol,
                               nodeAttr = gi.dc)


selected.str = genes%in%gi.str$ENSG
selectedGenes.str = genes[selected.str]
gi.str<-gi.str[match(selectedGenes.str,gi.str$ENSG),]

selectedAdj.str=adjacency.str[selected.str, selected.str];
dimnames(selectedAdj.str) = list(selectedGenes.str, selectedGenes.str)

allnodes<-gi[union(selectedGenes.str,selectedGenes.dc),]
cyt = exportNetworkToCytoscape(selectedAdj.str,
                               edgeFile = file.path(outdir,"Edges_ConsAdjStr.txt"),
                               nodeFile = file.path(outdir,"Nodes_ConsAdjStr.txt"),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = selectedGenes.str,
                               altNodeNames = gi.str$geneSymbol,
                               nodeAttr = gi.str)

cyt = exportNetworkToCytoscape(selectedAdj.str,
                               edgeFile = file.path(outdir,"Edges_ConsAdjStr.txt"),
                               nodeFile = file.path(outdir,"Nodes_All.txt"),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = union(selectedGenes.str,selectedGenes.dc),
                               altNodeNames = allnodes$geneSymbol,
                               nodeAttr = allnodes)
