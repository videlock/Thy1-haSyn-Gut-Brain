
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)

setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))


datlist<-readRDS("FinalProcessedData.rds")
an<-fData(readRDS("../1-RNAseqWorkflow/C_QCandNormalization/data/rawCountEset.rds"))
hasEGID<-rownames(an[!is.na(an$NCBI.gene.ID),])

grplist<-list(dc1m="Colon, 1 month",
              dc3m="Colon, 3 months",
              str1m="Striatum, 1 month",
              str3m="Striatum, 3 months"
)

# re-order samples by clustering within group
ReorderSamps<-function(datexpr,target){
  getInd<-function(subgrp){
    data1<-t(datexpr[,names(datexpr)%in%rownames(target[target$GT==subgrp,])])
    dend<-as.dendrogram(hclust(dist(data1)))
    Ind<-rev(order.dendrogram(dend))
    return(Ind)
  }
  samps1<-rownames(target[target$GT=="Hem",getInd("Hem")])
  samps2<-rownames(target[target$GT=="WT",getInd("WT")])
  samps<-c(samps1,samps2)
  datexpr<-datexpr[,samps]
  return(datexpr)
}


ann_colors = list(
  Genotype = c("Thy1-haSyn"="firebrick", WT="white"))

hmlist<-list()

for(grp in names(grplist)){
  datexpr<-ReorderSamps(datlist[[grp]]$dat.expr,datlist[[grp]]$target)
  target<-datlist[[grp]]$target[colnames(datexpr),]
  target$GT<-gsub("Hem","Thy1-haSyn",target$GT)
  tt<-datlist[[grp]]$tt[rownames(datlist[[grp]]$tt)%in%hasEGID,]
  hmlist[[grp]]<-pheatmap(mat = datexpr[rownames(tt)[1:50],],scale = "row",
                          cluster_cols = F,
                          annotation_col = data.frame(Genotype=target$GT, row.names = rownames(target)),
                          annotation_colors = ann_colors,
                          annotation_legend = F,
                          show_colnames = F,
                          labels_row = as.character(tt$Gene.name[1:50]),
                          breaks = seq(-2, 2, length.out = 101),
                          treeheight_row = 10,silent = F,
                          fontsize_row = 10,fontsize_col = 9,fontsize = 9,
                          main = grplist[[grp]],legend = F)
}



callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

hmlist2<-list()

for(grp in names(grplist)){
  datexpr<-ReorderSamps(datlist[[grp]]$dat.expr,datlist[[grp]]$target)
  target<-datlist[[grp]]$target[colnames(datexpr),]
  target$GT<-gsub("Hem","Thy1-haSyn",target$GT)
  tt<-datlist[[grp]]$tt[rownames(datlist[[grp]]$tt)%in%hasEGID,]
  hmlist2[[grp]]<-pheatmap(mat = datexpr[rownames(tt)[1:50],],scale = "row",
                          cluster_cols = F,clustering_callback = callback,
                          annotation_col = data.frame(Genotype=target$GT, row.names = rownames(target)),
                          annotation_colors = ann_colors,
                          annotation_legend = F,
                          show_colnames = F,
                          labels_row = as.character(tt$Gene.name[1:50]),
                          breaks = seq(-2, 2, length.out = 101),
                          treeheight_row = 10,silent = F,
                          fontsize_row = 8,fontsize_col = 8,fontsize = 7,
                          main = grplist[[grp]],legend = F)
}






AnnotLegendPlot<-pheatmap(mat = datexpr[rownames(tt)[1:50],],scale = "row",
                          cluster_cols = F,
                          annotation_col = data.frame(Genotype=target$GT, row.names = rownames(target)),
                          annotation_colors = ann_colors,
                          annotation_legend = T,
                          show_colnames = F,
                          labels_row = as.character(tt$Gene.name[1:50]),
                          breaks = seq(-2, 2, length.out = 101),
                          treeheight_row = 10,silent = F,
                          fontsize_row = 9,fontsize_col = 9,
                          main = grplist[[grp]],
                          legend = F)


MatLegendPlot<-pheatmap(mat = datexpr[rownames(tt)[1:50],],scale = "row",
                        cluster_cols = F,
                        annotation_col = data.frame(Genotype=target$GT, row.names = rownames(target)),
                        annotation_colors = ann_colors,annotation_names_row = F,
                        annotation_legend = F,annotation_names_col = F,
                        show_colnames = F,
                        labels_row = as.character(tt$Gene.name[1:50]),
                        breaks = seq(-2, 2, length.out = 101),
                        treeheight_row = 10,silent = F,
                        width = 2,height = 4,fontsize_row = 4,fontsize_col = 4,
                        main = grplist[[grp]],legend = T)




gl2=list(hmlist2$dc1m[[4]],hmlist2$dc3m[[4]],
        hmlist2$str1m[[4]],hmlist2$str3m[[4]],
        MatLegendPlot$gtable$grobs[[6]],
        AnnotLegendPlot$gtable$grobs[[7]])



pdf("GeneExpressionHeatmaps.pdf",height = 6.5*.9,width = 10*.9)
grid.arrange(
  grobs = gl2,
  widths = c(3, 3, 3, 3,1.6),
  heights=c(0.1,4,0.1,0.1,4),
  layout_matrix = rbind(c(1, 2, 3, 4, NA),
                        c(1, 2, 3, 4, 6),
                        c(1, 2, 3, 4, NA),
                        c(1, 2, 3, 4, NA),
                        c(1, 2, 3, 4, 5))
)
dev.off()
