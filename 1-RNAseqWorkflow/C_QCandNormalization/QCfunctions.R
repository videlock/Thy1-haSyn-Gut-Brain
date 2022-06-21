
# setup ---------
library(Biobase)
library(limma)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(MASS)
library(pheatmap)
library(tidyverse)
library(gridExtra)
library(grid)
library(WGCNA)


# setup directory structure ---------
dirSetup<-function(out.name){
  # outname is the name of the main directory
  # path.eset is the path for the rawCountEset set in B6_makeCountTable.R
  dir.create(out.name,showWarnings = FALSE)
  datdir<<-file.path(out.name,"data")
  figdir<<-file.path(out.name,"figs")
  dir.create(datdir,showWarnings = FALSE)
  dir.create(figdir,showWarnings = FALSE)
}


# load eset and restrict to group---------

loadEset<-function(Grp,path.eset,Outliers=NULL){
  eset<-readRDS(path.eset) 
  eset$Reads<-as.numeric(gsub(",","",eset$Reads))
  if(!is.null(Outliers)){
    eset<-eset[,!rownames(pData(eset))%in%Outliers]
  }
  if(Grp %in% c("dc","str")){
    eset<-eset[,eset$tissue==grp]
  }else{
    eset<-eset[,paste0(eset$tissue,eset$Time,"m")==grp]
  }
  return(eset)
}


# filter ---------------

filterEset<-function(Eset, CountsGrThan=5, inXpercent=30,
                     savePath=NULL){
  temp<-as.data.frame(exprs(Eset))
  temp$flag<-0
  for (i in 1:ncol(exprs(Eset))){
    temp$flag<-temp$flag+(temp[,i]>=CountsGrThan)
  }
  xPercent=round((inXpercent/100)*dim(Eset)[2],digits = 0)
  table(temp$flag>=xPercent) 
  keep=temp$flag>=xPercent
  filteredEset<-Eset[keep,]
  if(!is.null(savePath)){
    saveRDS(filteredEset,savePath)
  }
  return(filteredEset)
}


# outlier detection ----------

# for outlier detection a simple log2 transformation is made and rows are samples

getDatExprs<-function(Eset){
  dat.exprs<-exprs(Eset)
  dat.exprs2<-log2(dat.exprs+1)
  x2<-as.data.frame(t(dat.exprs2))
  rownames(x2)<-colnames(dat.exprs2)
  return(x2)
}

# WGCNA good genes and samples function
getGSG<-function(DatExprs){
  gsg<-goodSamplesGenes(DatExprs, minNGenes=1, verbose = 0)
  if(gsg$allOK){
    print("All samples and genes are good")
    return(DatExprs[gsg$goodSamples, gsg$goodGenes])
  }
}

# connectivity

findZKoutliers<-function(DatExprs, sdout=2, returnData=F){
  normadj <- (0.5+0.5*bicor(t(DatExprs), use='pairwise.complete.obs'))^2
  netsummary <- fundamentalNetworkConcepts(normadj); 
  K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
  C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
  outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
  print(paste(":There are ",sum(outliers),
              " outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep=""))
  outlierSamps<-colnames(t(DatExprs))[outliers]
  if(length(outlierSamps)>0){
    print(table(outliers))
    print(paste("Outliers:",outlierSamps)); 
  }
  
  if(returnData){
    keepSamples=!rownames(DatExprs)%in%outlierSamps
    return(DatExprs[keepSamples,])
  }else{
    return(outlierSamps)
  }
}

# sample tree
plotSampleTree<-function(DatExprs, main="Sample clustering to detect outliers",
                         savePath="SampleClustering.pdf"){
  sampleTree = hclust(dist(DatExprs), method = "average")
  pdf(savePath,height = 5,width = 8)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = main, sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  dev.off()
}


# density plots---------------
dlines<-function(counts,log=TRUE,main=""){
  if(log==TRUE){ 
    dat <- log(counts[,1],10)
    xlab="Raw read counts per gene (log10)"
  }else{
    dat<-counts[,1]
    xlab="Raw read counts per gene"
  }
  d <- density(dat)
  plot(d,xlim=c(1,8),main=main,ylim=c(0,.45),xlab=xlab, ylab="Density",cex.main=2)
  
  for (s in 2:dim(counts)[2]){
    if(log==TRUE){
      dat <- log(counts[,s],10) 
    }else{
      dat <- counts[,s]
    }
    d <- density(dat)
    lines(d)
  } 
}


dboxes<-function(counts, log=TRUE,main=""){
  if(log==TRUE){
    dat <- log(counts,10)
    ylab="Raw read counts per gene (log10)"
  }else{
    dat<-counts
    ylab="Raw read counts per gene"
  }
  suppressWarnings(boxplot(dat, main=main, xlab="", ylab=ylab,axes=FALSE))
  axis(2)
  axis(1,at=c(1:dim(counts)[2]),labels=colnames(counts),las=2,cex.axis=0.8)
  
}

printDensityAndMeanVarPlots<-function(fig.dir=".",
                                      rawData, FiltData, FiltNormData,
                                      DispData,
                                      normMethod="VSD",
                                      densityPlotName="densityplots.pdf",
                                      meanVarPlotName="meanvar.pdf"){
  pdf(file.path(fig.dir,densityPlotName),width = 8, height = 9)
  lm=rbind(c(7,1,2,2),
           c(8,3,4,4),
           c(9,5,6,6))
  layout(lm)
  dlines(rawData, main = paste0("n=",dim(rawData)[1]))
  dboxes(rawData)
  dlines(FiltData,main = paste0("n=",dim(FiltData)[1]))
  dboxes(FiltData)
  dlines(FiltNormData,log = F,main = paste0("n=",dim(FiltNormData)[1]))
  dboxes(FiltNormData,log = F)
  
  par(mar = c(0,0,0,0))
  plot.new();text(0.5,0.5,"raw",cex=3)
  plot.new();text(0.5,0.5,"filtered",cex=3)
  plot.new();text(0.5,0.5,normMethod,cex=3)
  dev.off()
  
  pdf(file.path(fig.dir,meanVarPlotName),width = 5, height = 5)
  plotDispEsts(DispData)
  dev.off()
}



# make target with picard variables -----------

makeTargetDat<-function(picard.path, Conds){
  picard<-readRDS(picard.path)
  picard<-picard[,c(
    "N_unmapped", "N_multimapping","N_noFeature","N_ambiguous", #star count info
    
    #AlignmentSummaryMetrics
    "PF_READS", "PF_HQ_ALIGNED_READS",
    "PF_INDEL_RATE",  "PF_HQ_ERROR_RATE", 
    "PF_ALIGNED_BASES", "PF_HQ_ALIGNED_BASES", "PF_HQ_ALIGNED_Q20_BASES", 
    "PF_MISMATCH_RATE", 
    
    #rnaseqmetrics
    "PF_BASES", "CODING_BASES", "PCT_CODING_BASES", 
    "UTR_BASES", "PCT_UTR_BASES",
    "INTRONIC_BASES", "PCT_INTRONIC_BASES",
    "INTERGENIC_BASES", "PCT_INTERGENIC_BASES", 
    "CORRECT_STRAND_READS", "PCT_CORRECT_STRAND_READS",
    "INCORRECT_STRAND_READS",
    "PCT_MRNA_BASES", "PCT_USABLE_BASES",
    "MEDIAN_CV_COVERAGE", 
    
    #DuplicationMetrics
    "UNPAIRED_READ_DUPLICATES",  "PERCENT_DUPLICATION")]
  
  picardgrp<-picard[rownames(Conds),]
  picardgrp<-as.data.frame(apply(picardgrp, 2, as.numeric))
  rownames(picardgrp)<-rownames(Conds)
  target<-cbind(Conds,picardgrp)
  n=dim(Conds)[2]
  target<-target[,!(colSums(is.na(target)) > 0)]
  temp1<-target[,-c(1:n)]
  temp1<-temp1[,colVars(as.matrix(temp1))!=0]
  temp1<-t(scale(temp1))
  temp1<-as.data.frame(t(temp1))
  target<-cbind(target[,1:n],temp1)
  return(target)
}




# add sequncing PCs to target--------
addSeqPCs<-function(Target, Conds, norm.Expr){
  normExpr<-as.data.frame(exprs(filtNormEset))
  if(unique(match(names(norm.Expr),Target$sample)==seq(1,nrow(Target),1))){
    tar1=Target[,c((1+dim(Conds)[2]):dim(Target)[2])]
    thisdat <- t(scale(tar1,scale=F))
    PC.metadata <- prcomp(thisdat,center=F);
    if(ncol(PC.metadata$rotation)>10){
      topPC1 <- PC.metadata$rotation[,1:10]
    }else{
      topPC1 <- PC.metadata$rotation[,1:ncol(PC.metadata$rotation)]
    }
    
    varexp <- (PC.metadata$sdev)^2 / sum(PC.metadata$sdev^2)
    seqtopvar <<- varexp[1:10]
    colnames(topPC1) <- paste("Seq.PC\n",colnames(topPC1)," (",signif(100*seqtopvar[1:7],2),"%)",sep="")
    Target$Seq.PC1=as.numeric(topPC1[,1])
    Target$Seq.PC2=as.numeric(topPC1[,2])
    Target$Seq.PC3=as.numeric(topPC1[,3])
    Target$Seq.PC4=as.numeric(topPC1[,4])
    Target$Seq.PC5=as.numeric(topPC1[,5])
    Target$Seq.PC6=as.numeric(topPC1[,6])
    Target$Seq.PC7=as.numeric(topPC1[,7])
    return(Target)
  }else{
    print("row names of target and column names of expression data must match")
  }
}

# expression PCs-----------------
PCvar<-function(data){
  thisdat.HTSC <- t(scale(t(data),scale=F))
  PC.HTSC <- prcomp(thisdat.HTSC,center=F);
  TopPC1 <- PC.HTSC$rotation[,1:5];
  varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
  topvar <- varexp[1:10]
  colnames(TopPC1) <- paste("Exp_",
                            colnames(TopPC1),"\n(", signif(100*topvar[1:5],2),"%)",sep="")
  return(list(TopPC1=TopPC1,
              varexp=varexp,
              topvar=topvar))
}


plotvar<-function(pclist,expPCs,
                  outdir=get("figdir", envir = parent.frame()),
                  fname){
  
  pdf(file.path(outdir,fname),height = 4,width = 7.5)
  par(mfrow=c(1,2))
  plot(pclist[[expPCs]][["topvar"]], xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", 
       ylab = "variance explained", main="Variance explained by Expression PCs",cex.main=1)
  plot(pclist$seq$topvar, xlim = c(0, 10), type = "b", pch = 16, xlab = "principal components", 
       ylab = "variance explained", main="Variance explained by sequencing PCs",cex.main=1)
  dev.off()
}

# MDS plots---------------------


getplotdata<-function(data,conds,kfit,topn=NULL){
  cv <- as.matrix(conds)
  dim(cv) <- c(1,prod(dim(cv)))
  names <- rownames(table(cv))
  col_default <- c("red","blue","green","orange","bisque4", "black",
                   "brown","cyan", "darkgreen", "darkgrey", "darkmagenta",
                   "darkolivegreen", "darkorange", "darkred", "darkslateblue",
                   "darkturquoise", "floralwhite", "greenyellow","grey", 
                   "lightcyan", "lightcyan1", "lightgreen", "lightsteelblue1",
                   "lightyellow","magenta", "mediumpurple3","midnightblue",
                   "paleturquoise", "pink", "plum1", "plum2", "royalblue",
                   "saddlebrown", "salmon", "sienna3", "skyblue", "skyblue3",
                   "steelblue", "tan", "thistle1", "thistle2", "turquoise",
                   "violet", "white", "yellowgreen", "grey60 ", "orangered4",
                   "brown4", "darkorange2", "ivory")
  col <- rainbow(length(names))
  names(col)<- names
  clab <- matrix(col[as.matrix(conds)],nrow=nrow(conds),ncol=ncol(conds))
  colnames(clab) <- colnames(conds) 
  rownames(clab)<-rownames(conds)
  n=dim(conds)[2]
  for(k in 1:n){
    a=1;for (i in unique(clab[,k])){
      clab[which(clab[,k]==i),k]=col_default[a]
      a=a+1
    }
  }
  if(!is.null(topn)){
    rowvar <- apply(data,1,var)
    ordV <- order(rowvar,decreasing=TRUE)
    topset <- head(ordV, n=topn)
    data <- data[topset,]
    main=paste0("Top ",topn, " by variance")
  }else{
    main="All genes"
  }
  ldat <- dist(t(data))
  invisible(fit <- isoMDS(ldat, k=kfit,))
  
  return(list(fit=fit,clab=clab, conds=conds, main=main))
}


makeplot<-function(plotdat,c1,c2,main=NULL,legpos="topleft",
                   legH=FALSE,legtxtw=0.045,legxint=0.25,
                   axes=T,frame=T,sizeText=get("sizeText",envir = parent.frame())){
  fit=plotdat$fit
  conds=plotdat$conds
  clab=plotdat$clab
  if(!is.null(main)){
    plot.title<-main
  }else{
    plot.title<-plotdat$main
  }
  x <- fit$points[,c1]; y <- fit$points[,c2]
  for(i in 2:dim(conds)[2]){
    plot(x, y, xlab=paste("Coordinate",c1),
         ylab=paste("Coordinate",c2),
         main=plot.title, type="p", pch=16, cex=1,
         col=clab[,i], xlim=c(min(x)*1.4, max(x)*1.4),
         ylim=c(min(y)*1.4, max(y)*1.4),axes=axes, frame.plot=frame)
    text(x,y+(max(y)-min(y))/50, labels = conds$sample, cex=sizeText,col=clab[,i])
    
    if(legH==TRUE){
      legend(legpos, legend=unique(as.factor(conds[,i])),fill=unique(clab[,i]),cex=sizeText, 
             horiz=TRUE,x.intersp=legxint, text.width=legtxtw,bty = "n")
    }else{
      legend(legpos, legend=unique(as.factor(conds[,i])),fill=unique(clab[,i]),cex=sizeText,bty = "n")
    }
    
  }
}


makemds.set<-function(data=get("normDat", envir = parent.frame()),
                      conds=get("conds", envir = parent.frame()),
                      kfit=6
){
  mdlist<-list(All = getplotdata(data = data, conds = conds, kfit = kfit,topn = NULL),
               top500 = getplotdata(data = data, conds = conds, kfit = kfit,topn = 500),
               top1000 = getplotdata(data = data, conds = conds, kfit = kfit,topn = 1000)
  )
  
  return(mdlist)
  
}

mdsplot1v2<-function(mds.set,
                     outdir=get("figdir", envir = parent.frame()),
                     fname){
  pdf(file.path(outdir,fname),width = 10, height = 5)
  sizeText=1; par(mfrow=c(1,3))
  for(i in 1:length(mds.set)){
    plotdat<-mds.set[[i]]
    makeplot(plotdat = plotdat,1,2,)
  }
  dev.off()
}

mdsplotAll<-function(mds.set,
                     outdir=get("figdir", envir = parent.frame()),
                     fname){
  pdf(file.path(outdir,fname),width = 10, height = 5)
  sizeText=.75; par(mar = c(0,0,0,0))
  layout(mat = matrix(seq(1,16),ncol = 4, byrow = T),
         heights = c(.5,4,4,4),widths = c(2,4,4,4))
  plot.new()
  plot.new();text(0.5,0.5,"All",cex=2)
  plot.new();text(0.5,0.5,"Top 500 (var)",cex=2)
  plot.new();text(0.5,0.5,"Top 1000 (var)",cex=2)
  plot.new();text(0.5,0.5,"1 v 2",cex=2)
  for(i in 1:length(mds.set)){
    plotdat<-mds.set[[i]]
    makeplot(plotdat = plotdat,1,2,main="",axes = F)
  }
  plot.new();text(0.5,0.5,"3 v 4",cex=2)
  for(i in 1:length(mds.set)){
    plotdat<-mds.set[[i]]
    makeplot(plotdat = plotdat,3,4,main="",axes = F)
  }
  plot.new();text(0.5,0.5,"5 v 6",cex=2)
  for(i in 1:length(mds.set)){
    plotdat<-mds.set[[i]]
    makeplot(plotdat = plotdat,5,6,main="",axes = F)
  }
  dev.off()
}


mdsCompPlot<-function(mds.list,
                      outdir=get("figdir", envir = parent.frame()),
                      fname, fw=10, fh=5){
  pdf(file.path(outdir,fname),width = fw, height = fh)
  sizeText=.75; par(mar = c(0,0,0,0))
  layout(mat = matrix(seq(1,(length(mds.list)+1)*4),ncol = 4,byrow = T),
         heights = c(.5,rep(4,length(mds.list))),widths = c(2,4,4,4))
  plot.new();text(0.5,0.5,"Coord 1 v 2",cex=2)
  plot.new();text(0.5,0.5,"All",cex=2)
  plot.new();text(0.5,0.5,"Top 500 (var)",cex=2)
  plot.new();text(0.5,0.5,"Top 1000 (var)",cex=2)
  
  for(i in 1:length(mds.list)){
    plot.new();text(0.5,0.5, names(mds.list)[i],cex=2)
    for(j in 1:length(mds.list[[i]])){
      plotdat<-mds.list[[i]][[j]]
      makeplot(plotdat = plotdat,1,2,main="",axes = F)
    }
  }
  dev.off()
}






# heatmaps----------------
get.r2mats<-function(topPC,Conds,Target){
  n=dim(Conds)[2]
  r2mat_meta = matrix(NA,nrow=5,ncol=n-1)
  r2mat_seq = matrix(NA,nrow=5,ncol=dim(Target)[2]-n)
  rownames(r2mat_meta) <- rownames(r2mat_seq) <- colnames(topPC)
  datMeta_model = Target[,2:n,drop=F]
  datSeq_model=Target[,-c(1:n)]
  
  colnames(r2mat_meta) <- colnames(datMeta_model)
  colnames(r2mat_seq) <- colnames(datSeq_model)
  for(i in c(1:5)){
    for(j in c(1:dim(datMeta_model)[2])){
      to_remove = which(is.na(datMeta_model[,j])==TRUE)
      if(length(to_remove)>0){
        tmp_topPC=topPC[-to_remove,i]
        tmp_meta=datMeta_model[-to_remove,j]}
      else{
        tmp_topPC=topPC[,i]
        tmp_meta=datMeta_model[,j]  
      }
      mod_mat=model.matrix(~tmp_meta)[,-1]
      mod=summary(lm(tmp_topPC~mod_mat))
      r2mat_meta[i,j]=mod$adj.r.squared
    }
  }   
  
  for(i in c(1:5)){
    for(j in c(1:dim(datSeq_model)[2])){
      to_remove = which(is.na(datSeq_model[,j])==TRUE)
      if(length(to_remove)>0){
        tmp_topPC=topPC[-to_remove,i]
        tmp_seq=datSeq_model[-to_remove,j]}
      else{
        tmp_topPC=topPC[,i]
        tmp_seq=datSeq_model[,j]  
      }
      mod=summary(lm(tmp_topPC~as.numeric(tmp_seq)))
      r2mat_seq[i,j]=mod$adj.r.squared
    }
  }
  
  return(list(r2mat_meta=r2mat_meta,
              r2mat_seq=r2mat_seq,
              datMeta_model=datMeta_model,
              datSeq_model=datSeq_model))
  
}


plothms<-function(pclist,expPCs,
                  Conds=get("conds", envir = parent.frame()),
                  Target=get("target", envir = parent.frame()), 
                  outdir=get("figdir", envir = parent.frame()),
                  fname){
  r2mats<-get.r2mats(pclist[[expPCs]][["TopPC1"]],Conds = Conds,Target = Target)
  hm1<-pheatmap(r2mats$r2mat_meta,cluster_cols = F,cluster_rows=FALSE, fontsize = 8, 
                silent=T,display_numbers = TRUE,fontsize_number = 8,cellwidth = 20, cellheight = 20)
  hm2<-pheatmap(t(r2mats$r2mat_seq),clustering_method="average",fontsize = 8,silent = T,
                cluster_cols=FALSE, display_numbers = TRUE,fontsize_number = 8, cellwidth = 40)
  hms<-arrangeGrob(grobs = list(hm1[[4]],hm2[[4]]), 
                   nrow = 1,widths = c(3,6))
  ggsave(plot = hms,filename = file.path(outdir,fname),
         width = 8,height = 5)
  
}

# correlations-----------


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y,use="pairwise.complete.obs",method="spearman"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.dots<-function(x,y,cex,col){
  points(x, y, pch = 19, col = cond, cex = 1)
  
}

panel.dotsSmall<-function(x,y,cex,col){
  points(x, y, pch = 19, col = cond, cex = 0.75)
  
}
  
# cex.points=1
cormatrix<-function(pclist,expPCs,
                    Pairsdat=get("pairsdat", envir = parent.frame()),
                    Cond=get("cond", envir = parent.frame()), 
                    outdir=get("figdir", envir = parent.frame()), fname, 
                    Narrow=F,
                    smallPoints=F){
  if(Narrow==TRUE){
    fw = 9; fh=8; cex.lab=1
    fname<-gsub("\\.pdf","N\\.pdf",fname)
  }else{
    fw = 15; fh=8; cex.lab=1.5
  }
  
  pdf(file.path(outdir,fname),width = fw, height=fh)
  
  if(smallPoints==T){
    pairs(cbind(pclist[[expPCs]][["TopPC1"]],Pairsdat),
          # col= cond,
          lower.panel = panel.dotsSmall,
          # pch=19,
          upper.panel = panel.cor,cex.labels=cex.lab,gap = .25)
  }else{
    pairs(cbind(pclist[[expPCs]][["TopPC1"]],Pairsdat),
          # col= cond,
          lower.panel = panel.dots,
          # pch=19,
          upper.panel = panel.cor,cex.labels=cex.lab,gap = .25)
  }
  
  dev.off()
}


correct<-function(modVars=get("modVars", envir = parent.frame()), 
                  Y=get("normExpr",envir = parent.frame()),
                  regColIdx=ifelse("Time"%in%modVars,4,3),interaction=FALSE){
  model<-paste("~",paste(modVars, collapse="+"))
  print(model)
  X =eval(parse(text = paste0("model.matrix(",model,")")))
  endIndex<-ifelse(interaction==FALSE,dim(X)[2],dim(X)[2]-1)
  print(head(X))
  print(paste0("cols to correct for: ",regColIdx,":",endIndex))
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  b = as.data.frame(t(beta))
  to_regress = (as.matrix(X[,regColIdx:endIndex,drop=FALSE]) %*% 
                  (as.matrix(beta[regColIdx:endIndex,,drop=FALSE])))
  datExprfit = normExpr - t(to_regress)
  return(datExprfit)
}



