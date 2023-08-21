
library(RRHO2)
library(limma)

setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))

outDir="rrhoFigs"
dir.create(outDir,showWarnings = F)


datlist<-readRDS("data/FinalProcessedData.rds")

datexpr.list<-list(Colon=datlist$dcAll$dat.expr,
                   Striatum=datlist$strAll$dat.expr)

f.list<-list(Colon=factor(datlist$dcAll$target$GTtime),
             Striatum=factor(datlist$strAll$target$GTtime))

an.list<-list(Colon=fData(readRDS("data/rawCountEset.rds"))[rownames(datexpr.list$Colon),],
              Striatum=fData(readRDS("data/rawCountEset.rds"))[rownames(datexpr.list$Striatum),]
)

tts<-list()

for(tis in c("Colon","Striatum")){
  GTtime<-f.list[[tis]]
  datexpr<-datexpr.list[[tis]]
  an<-an.list[[tis]]
  design<-model.matrix(~0+GTtime)
  colnames(design)<-levels(GTtime)
  contrasts<-makeContrasts(Hem_1 - WT_1, Hem_3 - WT_3,
                           Hem_3 - Hem_1, WT_3 - WT_1,
                           levels=design)
  fit <- lmFit(datexpr, design)
  fit$genes<-an[rownames(datexpr),c("Gene.name","Gene.description")]
  fit.cont<-contrasts.fit(fit, contrasts)
  fit.cont<-eBayes(fit.cont)
  tts[[paste0(tis,".GT.1m")]]<-topTable(fit.cont,coef = "Hem_1 - WT_1",
                                        adjust.method = "BH",number = "inf")
  tts[[paste0(tis,".GT.3m")]]<-topTable(fit.cont,coef = "Hem_3 - WT_3",
                                        adjust.method = "BH",number = "inf")
  tts[[paste0(tis,".Time.ASO")]]<-topTable(fit.cont,coef = "Hem_3 - Hem_1",
                                           adjust.method = "BH",number = "inf")
  tts[[paste0(tis,".Time.WT")]]<-topTable(fit.cont,coef = "WT_3 - WT_1",
                                          adjust.method = "BH",number = "inf")
  
}

# save(tts,file = "TopTables.rda")
# 
# load("TopTables.rda")


inputList<-lapply(tts, function(x){
    return(data.frame(gene = x$Gene.name,
                      rank = -1*log10(x$P.Value)*sign(x$t))
           )
})

# rm(tts)

# rrho - no adjustment --------------------
method = "hyper"
multipleTesting = "none"

rrhoList<-list()

rrhoList$GT.1m.CvS = 
  RRHO2_initialize(
    list1 = inputList$Colon.GT.1m[inputList$Colon.GT.1m$gene
                                  %in% 
                                    intersect(
                                      inputList$Colon.GT.1m$gene,
                                      inputList$Striatum.GT.1m$gene
                                    ),],
    list2 = inputList$Striatum.GT.1m[inputList$Striatum.GT.1m$gene
                                     %in% 
                                       intersect(
                                         inputList$Colon.GT.1m$gene,
                                         inputList$Striatum.GT.1m$gene
                                       ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (1m)", "Striatum (1m)")
  )


rrhoList$GT.3m.CvS = 
  RRHO2_initialize(
    list1 = inputList$Colon.GT.3m[inputList$Colon.GT.3m$gene
                                  %in% 
                                    intersect(
                                      inputList$Colon.GT.3m$gene,
                                      inputList$Striatum.GT.3m$gene
                                    ),],
    list2 = inputList$Striatum.GT.3m[inputList$Striatum.GT.3m$gene
                                     %in% 
                                       intersect(
                                         inputList$Colon.GT.3m$gene,
                                         inputList$Striatum.GT.3m$gene
                                       ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (3m)", "Striatum (3m)")
  )

rrhoList$GT.C.3mv1m = 
  RRHO2_initialize(
    list1 = inputList$Colon.GT.3m[inputList$Colon.GT.3m$gene
                                  %in% 
                                    intersect(
                                      inputList$Colon.GT.3m$gene,
                                      inputList$Colon.GT.1m$gene
                                    ),],
    list2 = inputList$Colon.GT.1m[inputList$Colon.GT.1m$gene
                                  %in% 
                                    intersect(
                                      inputList$Colon.GT.3m$gene,
                                      inputList$Colon.GT.1m$gene
                                    ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (3m)", "Colon (1m)")
  )


rrhoList$GT.S.3mv1m = 
  RRHO2_initialize(
    list1 = inputList$Striatum.GT.3m[inputList$Striatum.GT.3m$gene
                                     %in% 
                                       intersect(
                                         inputList$Striatum.GT.3m$gene,
                                         inputList$Striatum.GT.1m$gene
                                       ),],
    list2 = inputList$Striatum.GT.1m[inputList$Striatum.GT.1m$gene
                                     %in% 
                                       intersect(
                                         inputList$Striatum.GT.3m$gene,
                                         inputList$Striatum.GT.1m$gene
                                       ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Striatum (3m)", "Striatum (1m)")
  )

rrhoList$Time.ASO.CvS = 
  RRHO2_initialize(
    list1 = inputList$Colon.Time.ASO[inputList$Colon.Time.ASO$gene
                                     %in% 
                                       intersect(
                                         inputList$Colon.Time.ASO$gene,
                                         inputList$Striatum.Time.ASO$gene
                                       ),],
    list2 = inputList$Striatum.Time.ASO[inputList$Striatum.Time.ASO$gene
                                        %in% 
                                          intersect(
                                            inputList$Colon.Time.ASO$gene,
                                            inputList$Striatum.Time.ASO$gene
                                          ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (ASO)", "Striatum (ASO)")
  )

rrhoList$Time.WT.CvS = 
  RRHO2_initialize(
    list1 = inputList$Colon.Time.WT[inputList$Colon.Time.WT$gene
                                    %in% 
                                      intersect(
                                        inputList$Colon.Time.WT$gene,
                                        inputList$Striatum.Time.WT$gene
                                      ),],
    list2 = inputList$Striatum.Time.WT[inputList$Striatum.Time.WT$gene
                                       %in% 
                                         intersect(
                                           inputList$Colon.Time.WT$gene,
                                           inputList$Striatum.Time.WT$gene
                                         ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (WT)", "Striatum (WT)")
  )

rrhoList$Time.C.ASOvWT = 
  RRHO2_initialize(
    list1 = inputList$Colon.Time.ASO[inputList$Colon.Time.ASO$gene
                                     %in% 
                                       intersect(
                                         inputList$Colon.Time.ASO$gene,
                                         inputList$Colon.Time.WT$gene
                                       ),],
    list2 = inputList$Colon.Time.WT[inputList$Colon.Time.WT$gene
                                    %in% 
                                      intersect(
                                        inputList$Colon.Time.ASO$gene,
                                        inputList$Colon.Time.WT$gene
                                      ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (ASO)", "Colon (WT)")
  )

rrhoList$Time.S.ASOvWT = 
  RRHO2_initialize(
    list1 = inputList$Striatum.Time.ASO[inputList$Striatum.Time.ASO$gene
                                        %in% 
                                          intersect(
                                            inputList$Striatum.Time.ASO$gene,
                                            inputList$Striatum.Time.WT$gene
                                          ),],
    list2 = inputList$Striatum.Time.WT[inputList$Striatum.Time.WT$gene
                                       %in% 
                                         intersect(
                                           inputList$Striatum.Time.ASO$gene,
                                           inputList$Striatum.Time.WT$gene
                                         ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Striatum (ASO)", "Striatum (WT)")
  )


# rrho BY --------------
method = "hyper"
multipleTesting = "BY"


rrhoList.BY<-list()

rrhoList.BY$GT.1m.CvS = 
  RRHO2_initialize(
    list1 = inputList$Colon.GT.1m[inputList$Colon.GT.1m$gene
                                  %in% 
                                    intersect(
                                      inputList$Colon.GT.1m$gene,
                                      inputList$Striatum.GT.1m$gene
                                    ),],
    list2 = inputList$Striatum.GT.1m[inputList$Striatum.GT.1m$gene
                                     %in% 
                                       intersect(
                                         inputList$Colon.GT.1m$gene,
                                         inputList$Striatum.GT.1m$gene
                                       ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (1m)", "Striatum (1m)")
  )


rrhoList.BY$GT.3m.CvS = 
  RRHO2_initialize(
    list1 = inputList$Colon.GT.3m[inputList$Colon.GT.3m$gene
                                  %in% 
                                    intersect(
                                      inputList$Colon.GT.3m$gene,
                                      inputList$Striatum.GT.3m$gene
                                    ),],
    list2 = inputList$Striatum.GT.3m[inputList$Striatum.GT.3m$gene
                                     %in% 
                                       intersect(
                                         inputList$Colon.GT.3m$gene,
                                         inputList$Striatum.GT.3m$gene
                                       ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (3m)", "Striatum (3m)")
  )

rrhoList.BY$GT.C.3mv1m = 
  RRHO2_initialize(
    list1 = inputList$Colon.GT.3m[inputList$Colon.GT.3m$gene
                                  %in% 
                                    intersect(
                                      inputList$Colon.GT.3m$gene,
                                      inputList$Colon.GT.1m$gene
                                    ),],
    list2 = inputList$Colon.GT.1m[inputList$Colon.GT.1m$gene
                                  %in% 
                                    intersect(
                                      inputList$Colon.GT.3m$gene,
                                      inputList$Colon.GT.1m$gene
                                    ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (3m)", "Colon (1m)")
  )


rrhoList.BY$GT.S.3mv1m = 
  RRHO2_initialize(
    list1 = inputList$Striatum.GT.3m[inputList$Striatum.GT.3m$gene
                                     %in% 
                                       intersect(
                                         inputList$Striatum.GT.3m$gene,
                                         inputList$Striatum.GT.1m$gene
                                       ),],
    list2 = inputList$Striatum.GT.1m[inputList$Striatum.GT.1m$gene
                                     %in% 
                                       intersect(
                                         inputList$Striatum.GT.3m$gene,
                                         inputList$Striatum.GT.1m$gene
                                       ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Striatum (3m)", "Striatum (1m)")
  )

rrhoList.BY$Time.ASO.CvS = 
  RRHO2_initialize(
    list1 = inputList$Colon.Time.ASO[inputList$Colon.Time.ASO$gene
                                     %in% 
                                       intersect(
                                         inputList$Colon.Time.ASO$gene,
                                         inputList$Striatum.Time.ASO$gene
                                       ),],
    list2 = inputList$Striatum.Time.ASO[inputList$Striatum.Time.ASO$gene
                                        %in% 
                                          intersect(
                                            inputList$Colon.Time.ASO$gene,
                                            inputList$Striatum.Time.ASO$gene
                                          ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (ASO)", "Striatum (ASO)")
  )

rrhoList.BY$Time.WT.CvS = 
  RRHO2_initialize(
    list1 = inputList$Colon.Time.WT[inputList$Colon.Time.WT$gene
                                    %in% 
                                      intersect(
                                        inputList$Colon.Time.WT$gene,
                                        inputList$Striatum.Time.WT$gene
                                      ),],
    list2 = inputList$Striatum.Time.WT[inputList$Striatum.Time.WT$gene
                                       %in% 
                                         intersect(
                                           inputList$Colon.Time.WT$gene,
                                           inputList$Striatum.Time.WT$gene
                                         ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (WT)", "Striatum (WT)")
  )

rrhoList.BY$Time.C.ASOvWT = 
  RRHO2_initialize(
    list1 = inputList$Colon.Time.ASO[inputList$Colon.Time.ASO$gene
                                     %in% 
                                       intersect(
                                         inputList$Colon.Time.ASO$gene,
                                         inputList$Colon.Time.WT$gene
                                       ),],
    list2 = inputList$Colon.Time.WT[inputList$Colon.Time.WT$gene
                                    %in% 
                                      intersect(
                                        inputList$Colon.Time.ASO$gene,
                                        inputList$Colon.Time.WT$gene
                                      ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Colon (ASO)", "Colon (WT)")
  )

rrhoList.BY$Time.S.ASOvWT = 
  RRHO2_initialize(
    list1 = inputList$Striatum.Time.ASO[inputList$Striatum.Time.ASO$gene
                                        %in% 
                                          intersect(
                                            inputList$Striatum.Time.ASO$gene,
                                            inputList$Striatum.Time.WT$gene
                                          ),],
    list2 = inputList$Striatum.Time.WT[inputList$Striatum.Time.WT$gene
                                       %in% 
                                         intersect(
                                           inputList$Striatum.Time.ASO$gene,
                                           inputList$Striatum.Time.WT$gene
                                         ),],
    method = method, multipleTesting = multipleTesting,
    log10.ind=TRUE, boundary = 0.01,
    labels = c("Striatum (ASO)", "Striatum (WT)")
  )
# Warning message:
#   In RRHO2_initialize(list1 = inputList$Striatum.Time.ASO[inputList$Striatum.Time.ASO$gene %in%  :
#                                                             Inf was generated because of the multiple testing procedure. I.e., the multiple testing procedure cannot handle extreme small p-values. Suggest to use raw p-value (multipleTesting='none')






rrho2Plot<-function(obj1, obj2,
                    dir.text=c(rep(c("Upregulated","Downregulated"),2)),
                    plot.title,
                    file.name, outDir=".",
                    raster=FALSE,
                    res=300,
                    file.w=6.7, file.h=3.5,
                    cex.axistxt=0.8, fig.mar=c(3,3,0,2),
                    fig.oma=c(0,0,3,0),
                    widths=c(5.5, 5.5,1),
                    axtit.ml=1.5, dir.ml=0, dir.cex=0.5,
                    cex.title=1, ml.title=1.5,
                    mar.colorbar=c(3,3,0,0.5)
){
  color.bar <- function(lut, min, max=-min, 
                        nticks=11, 
                        ticks=seq(min, max, len=nticks), 
                        title='') {
    scale  <- (length(lut)-1)/(max-min)
    plot(c(0,10), c(min,max), type='n', bty='n', 
         xaxt='n', xlab='', yaxt='n', ylab='')
    mtext(title,2,2.3, cex=0.6)
    axis(2, round(ticks,0), las=1,cex.lab=0.5)
    for (i in 1:(length(lut)-1)) {
      y  <- (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
  }
  
  
  fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
    region <- match.arg(region, c("figure", "plot", "device"))
    pos <- match.arg(pos, c("topleft", "top", "topright", 
                            "left", "center", "right", 
                            "bottomleft", "bottom", "bottomright"))
    if(region %in% c("figure", "device")) {
      ds <- dev.size("in")
      # xy coordinates of device corners in user coordinates
      x <- grconvertX(c(0, ds[1]), from="in", to="user")
      y <- grconvertY(c(0, ds[2]), from="in", to="user")
      # fragment of the device we use to plot
      if(region == "figure") {
        # account for the fragment of the device that 
        # the figure is using
        fig <- par("fig")
        dx <- (x[2] - x[1])
        dy <- (y[2] - y[1])
        x <- x[1] + dx * fig[1:2]
        y <- y[1] + dy * fig[3:4]
      } 
    }
    # much simpler if in plotting region
    if(region == "plot") {
      u <- par("usr")
      x <- u[1:2]
      y <- u[3:4]
    }
    sw <- strwidth(text, cex=cex) * 60/100
    sh <- strheight(text, cex=cex) * 60/100
    x1 <- switch(pos,
                 topleft     =x[1] + sw, 
                 left        =x[1] + sw,
                 bottomleft  =x[1] + sw,
                 top         =(x[1] + x[2])/2,
                 center      =(x[1] + x[2])/2,
                 bottom      =(x[1] + x[2])/2,
                 topright    =x[2] - sw,
                 right       =x[2] - sw,
                 bottomright =x[2] - sw)
    y1 <- switch(pos,
                 topleft     =y[2] - sh,
                 top         =y[2] - sh,
                 topright    =y[2] - sh,
                 left        =(y[1] + y[2])/2,
                 center      =(y[1] + y[2])/2,
                 right       =(y[1] + y[2])/2,
                 bottomleft  =y[1] + sh,
                 bottom      =y[1] + sh,
                 bottomright =y[1] + sh)
    old.par <- par(xpd=NA)
    on.exit(par(old.par))
    text(x1, y1, text, cex=cex, ...)
    return(invisible(c(x,y)))
  }
  
  
  jet.colors  <- colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  hypermat.1 <- obj1$hypermat
  labels.1 <- obj1$labels
  method <- obj1$method
  maximum.1 <- max(hypermat.1,na.rm=TRUE)
  minimum.1 <- min(hypermat.1,na.rm=TRUE)
  hypermat.2 <- obj2$hypermat
  labels.2 <- obj2$labels
  maximum.2 <- max(hypermat.2,na.rm=TRUE)
  minimum.2 <- min(hypermat.2,na.rm=TRUE)
  maximum<-max(maximum.1,maximum.2)
  minimum<-max(minimum.1,minimum.2)
  
  if(substr(file.name,nchar(file.name)-3,nchar(file.name))==".png"){
    png(filename = file.path(outDir,file.name),width = file.w, 
        height = file.h, units = "in", res = res)
  }else{
    pdf(file = file.path(outDir,file.name),width = file.w, height = file.h)
  }
  
  colorGradient <- jet.colors(101)
  breaks <- seq(minimum,maximum,length.out = length(colorGradient) + 1)
  par(mar=fig.mar,oma=fig.oma)
  layoutmat=matrix(c(1,2,3),nrow = 1,ncol = 3,byrow = TRUE)
  layout(mat = layoutmat,widths = widths);
  image(hypermat.1, col = colorGradient,breaks=breaks,axes = FALSE,
        useRaster = raster)
  
  mtext(labels.1[2],2,axtit.ml,outer = F,cex = cex.axistxt)
  mtext(labels.1[1],1,axtit.ml,outer = F,cex = cex.axistxt)
  mtext(bquote(italic(.(dir.text[1]))),
        1,dir.ml,outer = F,cex = dir.cex, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text[1]))),
        3,dir.ml,outer = F,cex = dir.cex, 
        at=0.25,col = "red")
  
  mtext(bquote(italic(.(dir.text[2]))),
        1,0,outer = F,cex = dir.cex, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text[2]))),
        3,0,outer = F,cex = dir.cex, 
        at=0.75,col = "blue")
  
  mtext(bquote(italic(.(dir.text[3]))),
        2,0,outer = F,cex = dir.cex, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text[3]))),
        4,0,outer = F,cex = dir.cex, 
        at=0.25,col = "red")
  
  mtext(bquote(italic(.(dir.text[4]))),
        2,0,outer = F,cex = dir.cex, 
        at=0.75,col = "blue",)
  mtext(bquote(italic(.(dir.text[4]))),
        4,0,outer = F,cex = dir.cex, 
        at=0.75,col = "blue")
  
  fig_label("A.",pos = "bottomleft",cex = 1.5)
  image(hypermat.2, col = colorGradient,breaks=breaks,axes = FALSE,
        useRaster = raster)
  
  
  mtext(labels.2[2],2,axtit.ml,outer = F,cex = cex.axistxt)
  mtext(labels.2[1],1,axtit.ml,outer = F,cex = cex.axistxt)
  
  mtext(bquote(italic(.(dir.text[1]))),
        1,dir.ml,outer = F,cex = dir.cex, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text[1]))),
        3,dir.ml,outer = F,cex = dir.cex, 
        at=0.25,col = "red")
  
  mtext(bquote(italic(.(dir.text[2]))),
        1,0,outer = F,cex = dir.cex, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text[2]))),
        3,0,outer = F,cex = dir.cex, 
        at=0.75,col = "blue")
  
  mtext(bquote(italic(.(dir.text[3]))),
        2,0,outer = F,cex = dir.cex, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text[3]))),
        4,0,outer = F,cex = dir.cex, 
        at=0.25,col = "red")
  
  mtext(bquote(italic(.(dir.text[4]))),
        2,0,outer = F,cex = dir.cex, 
        at=0.75,col = "blue",)
  mtext(bquote(italic(.(dir.text[4]))),
        4,0,outer = F,cex = dir.cex, 
        at=0.75,col = "blue")

  fig_label("B.",pos = "bottomleft",cex = 1.5)
  
  mtext(plot.title,3,ml.title,outer = T,cex = cex.title)
  atitle ="-log10(P-value)"
  par(mar=mar.colorbar)
  color.bar(colorGradient, min = minimum, max = maximum, nticks = 7, 
            title = atitle)
  dev.off()
  
}
# aso v wt in colon v striatum---------------
obj1 = rrhoList$GT.1m.CvS
obj2 = rrhoList$GT.3m.CvS
dir.text = c("Upregulated - Colon", "Downregulated - Colon",
             "Upregulated - Striatum", "Downregulated - Striatum")
plot.title = "Overlap of Differential Expression (ASO vs WT) in Colon and Striatum"

file.name = "ColonVsStriatumGT.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)


obj1 = rrhoList.BY$GT.1m.CvS
obj2 = rrhoList.BY$GT.3m.CvS

file.name = "ColonVsStriatumGTBY.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)




# aso v wt in 1m v 3m---------------
obj1 = rrhoList$GT.C.3mv1m
obj2 = rrhoList$GT.S.3mv1m
dir.text = c("Upregulated - 1m", "Downregulated - 1m",
             "Upregulated - 3m", "Downregulated - 3m")
plot.title = "Overlap of Differential Expression (ASO vs WT) at 1 and 3 months"

file.name = "T1mVsT3mGT.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)



obj1 = rrhoList.BY$GT.C.3mv1m
obj2 = rrhoList.BY$GT.S.3mv1m

file.name = "T1mVsT3mGTBY.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)



# 1v3m in ASO by tissue ---------------

obj1 = rrhoList$Time.ASO.CvS
obj2 = rrhoList$Time.WT.CvS
dir.text = c("Upregulated - Colon", "Downregulated - Colon",
             "Upregulated - Striatum", "Downregulated - Striatum")
plot.title = "Overlap of Differential Expression (3m vs 1m) in Colon and Striatum"

file.name = "ColonVsStriatumTime.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)

file.name = "ColonVsStriatumTimeraster.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = T)

file.name = "ColonVsStriatumTime.png"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)

obj1 = rrhoList.BY$Time.ASO.CvS
obj2 = rrhoList.BY$Time.WT.CvS

file.name = "ColonVsStriatumTimeBY.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)

file.name = "ColonVsStriatumTimeBYraster.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = T)

file.name = "ColonVsStriatumTimeBY.png"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)



# 4 panel plot-------


rrho2Plotx4<-function(obj1, obj2,obj3,obj4,
                    dir.text.1=c(rep(c("Upregulated","Downregulated"),2)),
                    dir.text.2=c(rep(c("Upregulated","Downregulated"),2)),
                    set.title.1=NULL,set.title.2=NULL,
                    file.name, outDir=".",
                    labels.1=NULL,labels.2=NULL,
                    labels.3=NULL,labels.4=NULL,
                    plot.titles=NULL,
                    raster=FALSE,
                    res=300,
                    file.w=6.7, file.h=7,
                    fig.mar=c(3,3,3,2),
                    set.title.x=0.5,
                    set.title.y=0.25,
                    fig.oma=c(0,0,0,0),
                    widths=c(5.5, 5.5,1),
                    heights=c(1,5.5,1,5.5),
                    cex.dir=0.5,
                    cex.axistxt=0.8, 
                    cex.plot.title=1,
                    cex.set.title=1.5, 
                    ml.plot.title=1.5,
                    ml.axtit=1.5, ml.dir=0, 
                    mar.colorbar=c(3,3,0,0.5)
){
  color.bar <- function(lut, min, max=-min, 
                        nticks=11, 
                        ticks=seq(min, max, len=nticks), 
                        title='') {
    scale  <- (length(lut)-1)/(max-min)
    plot(c(0,10), c(min,max), type='n', bty='n', 
         xaxt='n', xlab='', yaxt='n', ylab='')
    mtext(title,2,2.3, cex=0.6)
    axis(2, round(ticks,0), las=1,cex.lab=0.5)
    for (i in 1:(length(lut)-1)) {
      y  <- (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
  }
  
  
  fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
    region <- match.arg(region, c("figure", "plot", "device"))
    pos <- match.arg(pos, c("topleft", "top", "topright", 
                            "left", "center", "right", 
                            "bottomleft", "bottom", "bottomright"))
    if(region %in% c("figure", "device")) {
      ds <- dev.size("in")
      # xy coordinates of device corners in user coordinates
      x <- grconvertX(c(0, ds[1]), from="in", to="user")
      y <- grconvertY(c(0, ds[2]), from="in", to="user")
      # fragment of the device we use to plot
      if(region == "figure") {
        # account for the fragment of the device that 
        # the figure is using
        fig <- par("fig")
        dx <- (x[2] - x[1])
        dy <- (y[2] - y[1])
        x <- x[1] + dx * fig[1:2]
        y <- y[1] + dy * fig[3:4]
      } 
    }
    # much simpler if in plotting region
    if(region == "plot") {
      u <- par("usr")
      x <- u[1:2]
      y <- u[3:4]
    }
    sw <- strwidth(text, cex=cex) * 60/100
    sh <- strheight(text, cex=cex) * 60/100
    x1 <- switch(pos,
                 topleft     =x[1] + sw, 
                 left        =x[1] + sw,
                 bottomleft  =x[1] + sw,
                 top         =(x[1] + x[2])/2,
                 center      =(x[1] + x[2])/2,
                 bottom      =(x[1] + x[2])/2,
                 topright    =x[2] - sw,
                 right       =x[2] - sw,
                 bottomright =x[2] - sw)
    y1 <- switch(pos,
                 topleft     =y[2] - sh,
                 top         =y[2] - sh,
                 topright    =y[2] - sh,
                 left        =(y[1] + y[2])/2,
                 center      =(y[1] + y[2])/2,
                 right       =(y[1] + y[2])/2,
                 bottomleft  =y[1] + sh,
                 bottom      =y[1] + sh,
                 bottomright =y[1] + sh)
    old.par <- par(xpd=NA)
    on.exit(par(old.par))
    text(x1, y1, text, cex=cex, ...)
    return(invisible(c(x,y)))
  }
  
  # setup ------
  hypermat.1 <- obj1$hypermat
  if(is.null(labels.1)){labels.1 <- obj1$labels}
  method <- obj1$method
  maximum.1 <- max(hypermat.1,na.rm=TRUE)
  minimum.1 <- min(hypermat.1,na.rm=TRUE)
  hypermat.2 <- obj2$hypermat
  if(is.null(labels.2)){labels.2 <- obj2$labels}
  maximum.2 <- max(hypermat.2,na.rm=TRUE)
  minimum.2 <- min(hypermat.2,na.rm=TRUE)
  maximum.set1<-max(maximum.1,maximum.2)
  minimum.set1<-max(minimum.1,minimum.2)
  hypermat.3 <- obj3$hypermat
  if(is.null(labels.3)){labels.3 <- obj3$labels}
  maximum.3 <- max(hypermat.3,na.rm=TRUE)
  minimum.3 <- min(hypermat.3,na.rm=TRUE)
  hypermat.4 <- obj4$hypermat
  if(is.null(labels.4)){labels.4 <- obj4$labels}
  maximum.4 <- max(hypermat.4,na.rm=TRUE)
  minimum.4 <- min(hypermat.4,na.rm=TRUE)
  maximum.set2<-max(maximum.3,maximum.4)
  minimum.set2<-max(minimum.3,minimum.4)
  jet.colors  <- colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", 
      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  colorGradient <- jet.colors(101)
  breaks.set1 <- seq(minimum.set1,maximum.set1,length.out = length(colorGradient) + 1)
  breaks.set2 <- seq(minimum.set2,maximum.set2,length.out = length(colorGradient) + 1)
  
  # layout -----
  
  if(substr(file.name,nchar(file.name)-3,nchar(file.name))==".png"){
    png(filename = file.path(outDir,file.name),width = file.w, 
        height = file.h, units = "in", res = res)
  }else{
    pdf(file = file.path(outDir,file.name),width = file.w, height = file.h)
  }
  
  par(mar=fig.mar,oma=fig.oma)
  layoutmat=matrix(c(1,1,0,2,3,4,5,5,0,6,7,8),nrow = 4,ncol = 3,byrow = TRUE)
  layout(mat = layoutmat,widths = widths, heights = heights)
  
  # title 1 -----
  
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(x=set.title.x,y=set.title.y,set.title.1,cex = cex.set.title,font=2)
  par(mar=fig.mar)
  # plot 1-------
  
  image(hypermat.1, col = colorGradient,breaks=breaks.set1,axes = FALSE,
        useRaster = raster)
  mtext(plot.titles[1],3,ml.plot.title,outer = F,cex = cex.plot.title)
  mtext(labels.1[2],2,ml.axtit,outer = F,cex = cex.axistxt)
  mtext(labels.1[1],1,ml.axtit,outer = F,cex = cex.axistxt)
  mtext(bquote(italic(.(dir.text.1[1]))),
        1,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.1[1]))),
        3,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.1[2]))),
        1,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text.1[2]))),
        3,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text.1[3]))),
        2,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.1[3]))),
        4,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.1[4]))),
        2,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue",)
  mtext(bquote(italic(.(dir.text.1[4]))),
        4,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  
  fig_label("A.",pos = "bottomleft",cex = 1.5)
  
  # plot 2----------
  image(hypermat.2, col = colorGradient,breaks=breaks.set1,axes = FALSE,
        useRaster = raster)
  mtext(plot.titles[2],3,ml.plot.title,outer = F,cex = cex.plot.title)
  mtext(labels.2[2],2,ml.axtit,outer = F,cex = cex.axistxt)
  mtext(labels.2[1],1,ml.axtit,outer = F,cex = cex.axistxt)
  
  mtext(bquote(italic(.(dir.text.1[1]))),
        1,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.1[1]))),
        3,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.1[2]))),
        1,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text.1[2]))),
        3,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text.1[3]))),
        2,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.1[3]))),
        4,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.1[4]))),
        2,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue",)
  mtext(bquote(italic(.(dir.text.1[4]))),
        4,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  
  fig_label("B.",pos = "bottomleft",cex = 1.5)
  
  # colorbar 1-----
  atitle ="-log10(P-value)"
  par(mar=mar.colorbar)
  color.bar(colorGradient, min = minimum.set1, max = maximum.set1, nticks = 7, 
            title = atitle)
  
  # title 2-----
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(x=set.title.x,y=set.title.y,set.title.2,cex = cex.set.title,font=2)
  par(mar=fig.mar)
  # plot 3-------
  image(hypermat.3, col = colorGradient,breaks=breaks.set2,axes = FALSE,
        useRaster = raster)
  mtext(plot.titles[3],3,ml.plot.title,outer = F,cex = cex.plot.title)
  mtext(labels.3[2],2,ml.axtit,outer = F,cex = cex.axistxt)
  mtext(labels.3[1],1,ml.axtit,outer = F,cex = cex.axistxt)
  mtext(bquote(italic(.(dir.text.2[1]))),
        1,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.2[1]))),
        3,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.2[2]))),
        1,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text.2[2]))),
        3,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text.2[3]))),
        2,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.2[3]))),
        4,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.2[4]))),
        2,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue",)
  mtext(bquote(italic(.(dir.text.2[4]))),
        4,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  
  fig_label("C.",pos = "bottomleft",cex = 1.5)
  
  # plot 4-------------
  image(hypermat.4, col = colorGradient,breaks=breaks.set2,axes = FALSE,
        useRaster = raster)
  mtext(plot.titles[4],3,ml.plot.title,outer = F,cex = cex.plot.title)
  mtext(labels.4[2],2,ml.axtit,outer = F,cex = cex.axistxt)
  mtext(labels.4[1],1,ml.axtit,outer = F,cex = cex.axistxt)
  
  mtext(bquote(italic(.(dir.text.2[1]))),
        1,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.2[1]))),
        3,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.2[2]))),
        1,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text.2[2]))),
        3,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  mtext(bquote(italic(.(dir.text.2[3]))),
        2,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.2[3]))),
        4,ml.dir,outer = F,cex = cex.dir, 
        at=0.25,col = "red")
  mtext(bquote(italic(.(dir.text.2[4]))),
        2,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue",)
  mtext(bquote(italic(.(dir.text.2[4]))),
        4,ml.dir,outer = F,cex = cex.dir, 
        at=0.75,col = "blue")
  
  fig_label("D.",pos = "bottomleft",cex = 1.5)
  
  # colorbar 2--------
  atitle ="-log10(P-value)"
  par(mar=mar.colorbar)
  color.bar(colorGradient, min = minimum.set2, max = maximum.set2, nticks = 7, 
            title = atitle)
  
  #-------
  dev.off()
  
}



obj1 = rrhoList$GT.1m.CvS
obj2 = rrhoList$GT.3m.CvS
obj3 = rrhoList$GT.C.3mv1m
obj4 = rrhoList$GT.S.3mv1m

dir.text.1 = c("Upregulated - Colon", "Downregulated - Colon",
             "Upregulated - Striatum", "Downregulated - Striatum")
dir.text.2 = c("Upregulated - 1m", "Downregulated - 1m",
             "Upregulated - 3m", "Downregulated - 3m")
set.title.1 = "Overlap in Colon and Striatum"
set.title.2 = "Overlap in 1 and 3 months"


file.name = "RRHO-unadjusted.png"
labels.1=c("Colon","Striatum")
labels.2=c("Colon","Striatum")
labels.3=c("1m","3m")
labels.4=c("1m","3m")
plot.titles=c("One month","Three months","Colon","Striatum")

rrho2Plotx4(obj1, obj2, obj3, obj4,
            dir.text.1,dir.text.2,
            set.title.1,set.title.2,
            file.name,
            outDir,
            labels.1 = labels.1,
            labels.2 = labels.2,
            labels.3 = labels.3,
            labels.4 = labels.4,
            plot.titles = plot.titles,
            raster = F)


obj1 = rrhoList.BY$GT.1m.CvS
obj2 = rrhoList.BY$GT.3m.CvS
obj3 = rrhoList.BY$GT.C.3mv1m
obj4 = rrhoList.BY$GT.S.3mv1m

file.name = "RRHO-BY.png"


rrho2Plotx4(obj1, obj2, obj3, obj4,
            dir.text.1,dir.text.2,
            set.title.1,set.title.2,
            file.name,
            outDir,
            labels.1 = labels.1,
            labels.2 = labels.2,
            labels.3 = labels.3,
            labels.4 = labels.4,
            plot.titles = plot.titles,
            raster = F)









# 1v3m in ASO by tissue ---------------

obj1 = rrhoList$Time.ASO.CvS
obj2 = rrhoList$Time.WT.CvS
dir.text = c("Upregulated - Colon", "Downregulated - Colon",
             "Upregulated - Striatum", "Downregulated - Striatum")
plot.title = "Overlap of Differential Expression (3m vs 1m) in Colon and Striatum"

file.name = "ColonVsStriatumTime.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)

file.name = "ColonVsStriatumTimeraster.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = T)

file.name = "ColonVsStriatumTime.png"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)

obj1 = rrhoList.BY$Time.ASO.CvS
obj2 = rrhoList.BY$Time.WT.CvS

file.name = "ColonVsStriatumTimeBY.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)

file.name = "ColonVsStriatumTimeBYraster.pdf"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = T)

file.name = "ColonVsStriatumTimeBY.png"
rrho2Plot(obj1, obj2, dir.text,plot.title,file.name,
          outDir,raster = F)


saveRDS(rrhoList,file = "rrhoListUnadjusted.rds")
saveRDS(rrhoList.BY,file = "rrhoListBY.rds")
saveRDS(inputList,file = "rrhoInputList.rds")
save(rrho2Plot,rrho2Plotx4, file = "rrho2Functions.rda")







# 1v3m by GT ---------------
plot.title="Overlap of Differential Expression (3m vs 1m) in Colon and Striatum"
file.name="ASOVsWTTime.pdf"
RRHO_obj.1<-rrhoList$Time.C.ASOvWT
RRHO_obj.2<-rrhoList$Time.S.ASOvWT
dirtxt.up.1="Upregulated - ASO"
dirtxt.down.1="Downregulated - ASO"
dirtxt.up.2="Upregulated - WT"
dirtxt.down.2="Downregulated - WT"
file.w=6.7
file.h=3.5
cex.axistxt=0.8
fig.mar=c(3,3,0,2)
fig.oma=c(0,0,3,0)
widths=c(5.5, 5.5,1)
axtit.ml=1.5
dir.ml=0
dir.cex=0.5
cex.title=1
ml.title=1.5
mar.colorbar=c(3,3,0,0.5)

hypermat.1 <- RRHO_obj.1$hypermat
labels.1 <- RRHO_obj.1$labels
method <- RRHO_obj.1$method
maximum.1 <- max(hypermat.1,na.rm=TRUE)
minimum.1 <- min(hypermat.1,na.rm=TRUE)

hypermat.2 <- RRHO_obj.2$hypermat
labels.2 <- RRHO_obj.2$labels
maximum.2 <- max(hypermat.2,na.rm=TRUE)
minimum.2 <- min(hypermat.2,na.rm=TRUE)
maximum<-max(maximum.1,maximum.2)
minimum<-max(minimum.1,minimum.2)


pdf(file = file.path(outDir,file.name),width = file.w, height = file.h)
colorGradient <- jet.colors(101)
breaks <- seq(minimum,maximum,length.out = length(colorGradient) + 1)
par(mar=fig.mar,oma=fig.oma)
layoutmat=matrix(c(1,2,3),nrow = 1,ncol = 3,byrow = TRUE)
layout(mat = layoutmat,widths = widths);
image(hypermat.1, col = colorGradient,breaks=breaks,axes = FALSE,)

mtext(labels.1[2],2,axtit.ml,outer = F,cex = cex.axistxt)
mtext(labels.1[1],1,axtit.ml,outer = F,cex = cex.axistxt)
mtext(bquote(italic(.(dirtxt.up.1))),
      1,dir.ml,outer = F,cex = dir.cex, 
      at=0.25,col = "red")
mtext(bquote(italic(.(dirtxt.up.1))),
      3,dir.ml,outer = F,cex = dir.cex, 
      at=0.25,col = "red")

mtext(bquote(italic(.(dirtxt.down.1))),
      1,0,outer = F,cex = dir.cex, 
      at=0.75,col = "blue")
mtext(bquote(italic(.(dirtxt.down.1))),
      3,0,outer = F,cex = dir.cex, 
      at=0.75,col = "blue")

mtext(bquote(italic(.(dirtxt.up.2))),
      2,0,outer = F,cex = dir.cex, 
      at=0.25,col = "red")
mtext(bquote(italic(.(dirtxt.up.2))),
      4,0,outer = F,cex = dir.cex, 
      at=0.25,col = "red")

mtext(bquote(italic(.(dirtxt.down.2))),
      2,0,outer = F,cex = dir.cex, 
      at=0.75,col = "blue",)
mtext(bquote(italic(.(dirtxt.down.2))),
      4,0,outer = F,cex = dir.cex, 
      at=0.75,col = "blue")

fig_label("A.",pos = "bottomleft",cex = 1.5)
image(hypermat.2, col = colorGradient,breaks=breaks,axes = FALSE)


mtext(labels.2[2],2,axtit.ml,outer = F,cex = cex.axistxt)
mtext(labels.2[1],1,axtit.ml,outer = F,cex = cex.axistxt)

mtext(bquote(italic(.(dirtxt.up.1))),
      1,dir.ml,outer = F,cex = dir.cex, 
      at=0.25,col = "red")
mtext(bquote(italic(.(dirtxt.up.1))),
      3,dir.ml,outer = F,cex = dir.cex, 
      at=0.25,col = "red")

mtext(bquote(italic(.(dirtxt.down.1))),
      1,0,outer = F,cex = dir.cex, 
      at=0.75,col = "blue")
mtext(bquote(italic(.(dirtxt.down.1))),
      3,0,outer = F,cex = dir.cex, 
      at=0.75,col = "blue")

mtext(bquote(italic(.(dirtxt.up.2))),
      2,0,outer = F,cex = dir.cex, 
      at=0.25,col = "red")
mtext(bquote(italic(.(dirtxt.up.2))),
      4,0,outer = F,cex = dir.cex, 
      at=0.25,col = "red")

mtext(bquote(italic(.(dirtxt.down.2))),
      2,0,outer = F,cex = dir.cex, 
      at=0.75,col = "blue",)
mtext(bquote(italic(.(dirtxt.down.2))),
      4,0,outer = F,cex = dir.cex, 
      at=0.75,col = "blue")

fig_label("B.",pos = "bottomleft",cex = 1.5)

mtext(plot.title,3,ml.title,outer = T,cex = cex.title)
atitle ="-log10(P-value)"
par(mar=mar.colorbar)
color.bar(colorGradient, min = minimum, max = maximum, nticks = 7, 
          title = atitle)
dev.off()


