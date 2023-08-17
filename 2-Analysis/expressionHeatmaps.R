
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)

setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))


datlist<-readRDS("data/FinalProcessedData.rds")
an<-fData(readRDS("data/rawCountEset.rds"))

hasEGID<-rownames(an[!is.na(an$NCBI.gene.ID),])

grplist<-list(dc1m="Gut, 1 month",
              dc3m="Gut, 3 months",
              str1m="Brain, 1 month",
              str3m="Brain, 3 months"
)

# re-order samples by clustering within group
ReorderSamps<-function(datexpr,target){
  getInd<-function(subgrp){
    data1<-t(datexpr[,names(datexpr)%in%rownames(target[target$GT==subgrp,])])
    dend<-as.dendrogram(hclust(dist(data1)))
    # dend <- reorder(dend, rowMeans(data1))
    Ind<-rev(order.dendrogram(dend))
    return(Ind)
  }
  samps1<-rownames(target[target$GT=="Hem",getInd("Hem")])
  samps2<-rownames(target[target$GT=="WT",getInd("WT")])
  samps<-c(samps1,samps2)
  datexpr<-datexpr[,samps]
  return(datexpr)
}




fontsize_row=6
fontsize_col=6
fontsize=7


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
                          annotation_col = data.frame(Genotype=target$GT, 
                                                      row.names = rownames(target) ),
                          annotation_colors = ann_colors,
                          annotation_legend = F,
                          show_colnames = F,
                          labels_row = as.character(tt$Gene.name[1:50]),
                          breaks = seq(-2, 2, length.out = 101),
                          treeheight_row = 10,silent = F,
                          fontsize_row = fontsize_row,fontsize_col = fontsize_col,fontsize = fontsize,
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
                           annotation_col = data.frame(Genotype=target$GT, 
                                                       row.names = rownames(target) ),
                           annotation_colors = ann_colors,
                           annotation_legend = F,
                           show_colnames = F,
                           labels_row = as.character(tt$Gene.name[1:50]),
                           breaks = seq(-2, 2, length.out = 101),
                           treeheight_row = 10,silent = F,
                           fontsize_row = fontsize_row,fontsize_col = fontsize_col,fontsize = fontsize,
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
                          fontsize_row = fontsize_row,fontsize_col = fontsize_col,fontsize = fontsize,
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
                        width = 2,height = 4,fontsize_row = fontsize_row,fontsize_col = fontsize_col,fontsize = fontsize,
                        main = grplist[[grp]],legend = T)




gl2=list(hmlist2$dc1m[[4]],hmlist2$dc3m[[4]],
         hmlist2$str1m[[4]],hmlist2$str3m[[4]],
         MatLegendPlot$gtable$grobs[[6]],
         AnnotLegendPlot$gtable$grobs[[7]])


a.ratio=10/6.5

w=6.5
h=w/a.ratio

pdf("expressionHeatmap.pdf",height = h,width = w)
grid.arrange(
  grobs = gl2,
  widths = c(3, 3, 3, 3,1.6),
  heights=c(0.1,2,0.1,1,1),
  layout_matrix = rbind(c(1, 2, 3, 4, NA),
                        c(1, 2, 3, 4, 6),
                        c(1, 2, 3, 4, NA),
                        c(1, 2, 3, 4, 5),
                        c(1, 2, 3, 4, NA))
)
dev.off()
