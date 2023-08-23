
load("regInputData.rda")

source("myMultivarFunctions.R")
library(sjPlot)
library(easystats)
library(report)
library(limma)
library(MASS)
library(ggplot2)
library(corrplot)
library(knitr)
library(broom)
library(flextable)

# Striatum --------

mod.age<-NULL

for(mod in names(str.MEs)){
  s.aso<-summary(lm(str.dat[,mod]~str.dat$Age))
  
  mod.age<-rbind(mod.age,
                 s.aso$coefficients[2,c("Estimate","Pr(>|t|)")])
  
}


mod.age<-as.data.frame(mod.age)
rownames(mod.age)<-names(str.MEs)
agemods<-rownames(mod.age[mod.age$`Pr(>|t|)`<0.05,])
agemods.str<-agemods
univ.fit<-NULL

for(mod in names(str.MEs)){
  
  if(mod %in% agemods){
    s.aso<-summary(lm(str.dat$ASO~str.dat[,mod]+str.dat$Age+str.dat[,mod]*str.dat$Age))
  }else{
    s.aso<-summary(lm(str.dat$ASO~str.dat[,mod]))
    
  }
  
  univ.fit<-rbind(univ.fit,
                  s.aso$coefficients[2,c("Estimate","Pr(>|t|)")])
  
}


univ.fit<-as.data.frame(univ.fit)
rownames(univ.fit)<-names(str.MEs)

uvmods.str<-rownames(univ.fit[univ.fit$`Pr(>|t|)`<0.05,])

str.uv.tab<-data.frame(Predictors=gsub("str\\.","",rownames(univ.fit[uvmods.str,])),
                       Estimates=univ.fit[uvmods.str,"Estimate"],
                       p=univ.fit[uvmods.str,"Pr(>|t|)"],row.names = uvmods.str)

for(i in 1:nrow(str.uv.tab)){
  if(rownames(str.uv.tab)[i]%in%agemods){
    str.uv.tab$Predictors[i]<-paste0(str.uv.tab$Predictors[i],"*")
  }
}

supMatTabs<-list(dfs=list(uv.str=str.uv.tab),
                 titles=list(uv.str="Striatum univariate predictors of Thy1-haSyn status"))

cor1 = cor(str.dat[,uvmods.str])
corrplot.mixed(cor1, lower.col = "black", number.cex = .7,
               upper.col=topo.colors(50))

if(length(intersect(uvmods.str,agemods.str))>0){
  fit.mv<-eval(parse(
    text = paste("lm(ASO ~", 
                 paste(uvmods.str,collapse = " + "),
                 "+Age+",
                 paste(paste(uvmods.str[uvmods.str%in%agemods],"Age",sep = "*"),collapse = " + "),
                 
                 ", data = str.dat)"),
  ))
}else{
  fit.mv<-eval(parse(
    text = paste("lm(ASO ~", 
                 paste(uvmods.str,collapse = " + "),
                 "+Age",
                 
                 ", data = str.dat)"),
  ))
}

sum.mv<-summary(fit.mv)
fit.mv.step<-stepAIC(fit.mv,direction = "backward",trace = FALSE)

restab<-as.data.frame(summary(fit.mv.step)$coefficients)
mvmods.str<-rownames(restab[restab$`Pr(>|t|)`<0.05,])
mvmods.str<-mvmods.str[!mvmods.str%in%c("(Intercept)","Age")]

str.mv.tab<-data.frame(Predictors=gsub("str\\.","",rownames(restab[mvmods.str,])),
                       Estimates=restab[mvmods.str,"Estimate"],
                       p=restab[mvmods.str,"Pr(>|t|)"],row.names = mvmods.str)

mvfit.str<-fit.mv.step

# Colon----------

mod.age<-NULL

for(mod in names(dc.MEs)){
  s.aso<-summary(lm(dc.dat[,mod]~dc.dat$Age))
  
  mod.age<-rbind(mod.age,
                 s.aso$coefficients[2,c("Estimate","Pr(>|t|)")])
  
}

mod.age<-as.data.frame(mod.age)
rownames(mod.age)<-names(dc.MEs)
agemods<-rownames(mod.age[mod.age$`Pr(>|t|)`<0.05,])
agemods.dc<-agemods

univ.fit<-NULL

for(mod in names(dc.MEs)){
  
  if(mod %in% agemods){
    s.aso<-summary(lm(dc.dat$ASO~dc.dat[,mod]+dc.dat$Age+dc.dat[,mod]*dc.dat$Age))
  }else{
    s.aso<-summary(lm(dc.dat$ASO~dc.dat[,mod]))
    
  }
  
  s.aso1m<-summary(lm(dc.dat$ASO[dc.dat$Time==1]~dc.dat[dc.dat$Time==1,mod]))
  s.aso3m<-summary(lm(dc.dat$ASO[dc.dat$Time==3]~dc.dat[dc.dat$Time==3,mod]))
  
  univ.fit<-rbind(univ.fit,
                  c(s.aso$coefficients[2,c("Estimate","Pr(>|t|)")],
                    s.aso1m$coefficients[2,c("Estimate","Pr(>|t|)")],
                    s.aso3m$coefficients[2,c("Estimate","Pr(>|t|)")]
                  ))
}



univ.fit<-as.data.frame(univ.fit)
names(univ.fit)<-c("asoEst","asoP","aso1mEst","aso1mP","aso3mEst","aso3mP")
rownames(univ.fit)<-names(dc.MEs)
uvmods.dc<-rownames(univ.fit[univ.fit$asoP<0.05,])
uvmods.dc1m<-rownames(univ.fit[univ.fit$aso1mP<0.05,])
uvmods.dc3m<-rownames(univ.fit[univ.fit$aso3mP<0.05,])

cor1 = cor(dc.dat[dc.dat$Time==1,uvmods.dc])
corrplot.mixed(cor1, lower.col = "black", number.cex = .7,
               upper.col=topo.colors(50))

dc.uv.tab<-data.frame(Predictors=gsub("dc\\.","",rownames(univ.fit[uvmods.dc,])),
                      Estimates=univ.fit[uvmods.dc,"asoEst"],
                      p=univ.fit[uvmods.dc,"asoP"],row.names = uvmods.dc)

for(i in 1:nrow(dc.uv.tab)){
  if(rownames(dc.uv.tab)[i]%in%agemods){
    dc.uv.tab$Predictors[i]<-paste0(dc.uv.tab$Predictors[i],"*")
  }
}


supMatTabs$dfs$uv.dc<-dc.uv.tab
supMatTabs$titles$uv.dc<-"Colon univariate predictors of Thy1-haSyn status"

supMatTabs$dfs$uv.dc1m<-data.frame(Predictors=gsub("dc\\.","",rownames(univ.fit[uvmods.dc1m,])),
                                   Estimates=univ.fit[uvmods.dc1m,"aso1mEst"],
                                   p=univ.fit[uvmods.dc1m,"aso1mP"],row.names = uvmods.dc1m)
supMatTabs$titles$uv.dc1m<-"Colon univariate predictors of Thy1-haSyn status (one month)"

supMatTabs$dfs$uv.dc3m<-data.frame(Predictors=gsub("dc\\.","",rownames(univ.fit[uvmods.dc3m,])),
                                   Estimates=univ.fit[uvmods.dc3m,"aso3mEst"],
                                   p=univ.fit[uvmods.dc3m,"aso3mP"],row.names = uvmods.dc3m)
supMatTabs$titles$uv.dc3m<-"Colon univariate predictors of Thy1-haSyn status (three months)"


uvmods.dc.2<-uvmods.dc[!uvmods.dc%in%c("dc.darkgrey","dc.cyan")]


fit.mv<-eval(parse(
  text = paste("lm(ASO ~",
               paste(uvmods.dc.2,collapse = " + "),
               "+Age+",
               paste(
                 paste(uvmods.dc.2[uvmods.dc.2%in%agemods],"Age",sep = "*"),
                 collapse = " + "),
               ", data = dc.dat)")
))

fit.mv.step<-stepAIC(fit.mv,direction = "backward",trace = FALSE)
sum.mv.step<-summary(fit.mv.step)
mvfit.dc<-fit.mv.step


# colon 1 month

cor1 = cor(dc.dat[dc.dat$Time==1,uvmods.dc1m])
corrplot.mixed(cor1, lower.col = "black", number.cex = .7,
               upper.col=topo.colors(50))

uvmods.dc1m.2<-uvmods.dc1m[!uvmods.dc1m%in%c("dc.darkgrey","dc.cyan")]


fit.mv.1m<-eval(parse(
  text = paste("lm(",
               paste("ASO ~",
                     paste(uvmods.dc1m.2,collapse = " + ")),
               ", data = dc.dat[dc.dat$Time==1,])")))

fit.mv.1m.step<-stepAIC(fit.mv.1m,direction = "backward",trace = T)
sum.mv.1m.step<-summary(fit.mv.1m.step)

restab<-as.data.frame(sum.mv.1m.step$coefficients)
mvmods.dc1m<-rownames(restab[restab$`Pr(>|t|)`<0.05,])
mvmods.dc1m<-mvmods.dc1m[!mvmods.dc1m%in%"(Intercept)"]

mvfit.dc1m<-fit.mv.1m.step


# colon 3 months
fit.mv.3m<-eval(parse(
  text = paste("lm(",
               paste("ASO ~",
                     paste(uvmods.dc3m,collapse = " + ")),
               ", data = dc.dat[dc.dat$Time==3,])")))
sum.mv.3m<-summary(fit.mv.3m)
mvfit.dc3m<-fit.mv.3m

restab<-as.data.frame(sum.mv.3m$coefficients)
mvmods.dc3m<-rownames(restab[restab$`Pr(>|t|)`<0.05,])
mvmods.dc3m<-mvmods.dc3m[!mvmods.dc3m%in%"(Intercept)"]


# Combined striatum and colon------------

# striatum matched samples
mod.age<-NULL

for(mod in names(str.MEs)){
  s.aso<-summary(lm(bg.dat[,mod]~bg.dat$Age))
  
  mod.age<-rbind(mod.age,
                 s.aso$coefficients[2,c("Estimate","Pr(>|t|)")])
  
}


mod.age<-as.data.frame(mod.age)
rownames(mod.age)<-names(str.MEs)
agemods<-rownames(mod.age[mod.age$`Pr(>|t|)`<0.05,])
agemods.str.bg<-agemods

univ.fit.bg<-NULL

for(mod in names(str.MEs)){
  
  if(mod %in% agemods){
    s.aso<-summary(lm(bg.dat$ASO~bg.dat[,mod]+bg.dat$Age+bg.dat[,mod]*bg.dat$Age))
  }else{
    s.aso<-summary(lm(bg.dat$ASO~bg.dat[,mod]))
    
  }
  
  univ.fit.bg<-rbind(univ.fit.bg,
                     s.aso$coefficients[2,c("Estimate","Pr(>|t|)")])
  
}

univ.fit.bg<-as.data.frame(univ.fit.bg)
rownames(univ.fit.bg)<-names(str.MEs)

uvmods.str.bg<-rownames(univ.fit.bg[univ.fit.bg$`Pr(>|t|)`<0.05,])

uv.str.bg<-data.frame(Predictors=gsub("str\\.","",rownames(univ.fit.bg[uvmods.str.bg,])),
                      Estimates=univ.fit.bg[uvmods.str.bg,"Estimate"],
                      p=univ.fit.bg[uvmods.str.bg,"Pr(>|t|)"],row.names = uvmods.str.bg)

for(i in 1:nrow(uv.str.bg)){
  if(rownames(uv.str.bg)[i]%in%agemods){
    uv.str.bg$Predictors[i]<-paste0(uv.str.bg$Predictors[i],"*")
  }
}


supMatTabs$dfs$uv.str.bg<-uv.str.bg
supMatTabs$titles$uv.str.bg<-"Striatum (matched samples) univariate predictors of Thy1-haSyn status"

if(length(uvmods.str.bg[uvmods.str.bg%in%agemods])>0){
  fit.mv.bg<-eval(parse(
    text = paste("lm(ASO ~",
                 paste(uvmods.str.bg,collapse = " + "),
                 "+Age+",
                 paste(
                   paste(uvmods.str.bg[uvmods.str.bg%in%agemods],"Age",sep = "*"),
                   collapse = " + "),
                 ", data = bg.dat)")
  ))
}else{
  fit.mv.bg<-eval(parse(
    text = paste("lm(ASO ~",
                 paste(uvmods.str.bg,collapse = " + "),
                 ", data = bg.dat)")
  ))
}

# summary(fit.mv.bg)
fit.mv.bg.step<-stepAIC(fit.mv.bg,direction = "backward",trace = FALSE)
mvfit.str.bg<-fit.mv.bg.step

restab<-as.data.frame(summary(fit.mv.bg.step)$coefficients)
mvmods.str.bg<-rownames(restab[restab$`Pr(>|t|)`<0.05,])
mvmods.str.bg<-mvmods.str.bg[!mvmods.str.bg%in%"(Intercept)"]


cor1 = cor(bg.dat[,uvmods.str.bg])
corrplot.mixed(cor1, lower.col = "black", number.cex = .7,
               upper.col=topo.colors(50))

# colon matched samples

mod.age<-NULL

for(mod in names(dc.MEs)){
  s.aso<-summary(lm(bg.dat[,mod]~bg.dat$Age))
  
  mod.age<-rbind(mod.age,
                 s.aso$coefficients[2,c("Estimate","Pr(>|t|)")])
  
}

mod.age<-as.data.frame(mod.age)
rownames(mod.age)<-names(dc.MEs)
agemods<-rownames(mod.age[mod.age$`Pr(>|t|)`<0.05,])
agemods.dc.bg<-agemods

univ.fit.bg<-NULL

for(mod in names(dc.MEs)){
  
  if(mod %in% agemods){
    s.aso<-summary(lm(bg.dat$ASO~bg.dat[,mod]+bg.dat$Age+bg.dat[,mod]*bg.dat$Age))
  }else{
    s.aso<-summary(lm(bg.dat$ASO~bg.dat[,mod]))
    
  }
  
  s.aso1m<-summary(lm(bg.dat$ASO[bg.dat$Time==1]~bg.dat[bg.dat$Time==1,mod]))
  s.aso3m<-summary(lm(bg.dat$ASO[bg.dat$Time==3]~bg.dat[bg.dat$Time==3,mod]))
  
  univ.fit.bg<-rbind(univ.fit.bg,
                     c(s.aso$coefficients[2,c("Estimate","Pr(>|t|)")],
                       s.aso1m$coefficients[2,c("Estimate","Pr(>|t|)")],
                       s.aso3m$coefficients[2,c("Estimate","Pr(>|t|)")]
                     ))
}


univ.fit.bg<-as.data.frame(univ.fit.bg)
names(univ.fit.bg)<-c("asoEst","asoP","aso1mEst","aso1mP","aso3mEst","aso3mP")
rownames(univ.fit.bg)<-names(dc.MEs)

uvmods.dc.bg<-rownames(univ.fit.bg[univ.fit.bg$asoP<0.05,])
uvmods.dc1m.bg<-rownames(univ.fit.bg[univ.fit.bg$aso1mP<0.05,])
uvmods.dc3m.bg<-rownames(univ.fit.bg[univ.fit.bg$aso3mP<0.05,])

cor1 = cor(bg.dat[,uvmods.dc.bg])
corrplot.mixed(cor1, lower.col = "black", number.cex = .7,
               upper.col=topo.colors(50))

uv.dc.bg<-data.frame(Predictors=gsub("dc\\.","",rownames(univ.fit.bg[uvmods.dc.bg,])),
                     Estimates=univ.fit.bg[uvmods.dc.bg,"asoEst"],
                     p=univ.fit.bg[uvmods.dc.bg,"asoP"],row.names = uvmods.dc.bg)

for(i in 1:nrow(uv.dc.bg)){
  if(rownames(uv.dc.bg)[i]%in%agemods){
    uv.dc.bg$Predictors[i]<-paste0(uv.dc.bg$Predictors[i],"*")
  }
}

supMatTabs$dfs$uv.dc.bg<-uv.dc.bg
supMatTabs$titles$uv.dc.bg<-"Colon (matched samples) univariate predictors of Thy1-haSyn status"

supMatTabs$dfs$uv.dc1m.bg<-data.frame(Predictors=gsub("dc\\.","",rownames(univ.fit.bg[uvmods.dc1m.bg,])),
                                      Estimates=univ.fit.bg[uvmods.dc1m.bg,"aso1mEst"],
                                      p=univ.fit.bg[uvmods.dc1m.bg,"aso1mP"],row.names = uvmods.dc1m.bg)
supMatTabs$titles$uv.dc1m.bg<-"Colon (matched samples) univariate predictors of Thy1-haSyn status (one month)"

supMatTabs$dfs$uv.dc3m.bg<-data.frame(Predictors=gsub("dc\\.","",rownames(univ.fit.bg[uvmods.dc3m.bg,])),
                                      Estimates=univ.fit.bg[uvmods.dc3m.bg,"aso3mEst"],
                                      p=univ.fit.bg[uvmods.dc3m.bg,"aso3mP"],row.names = uvmods.dc3m.bg)
supMatTabs$titles$uv.dc3m.bg<-"Colon (matched samples) univariate predictors of Thy1-haSyn status (three months)"

uvmods.dc1m.bg.2<-uvmods.dc1m.bg[!uvmods.dc1m.bg%in%c("dc.cyan")]

if(length(uvmods.dc1m.bg.2[uvmods.dc1m.bg.2%in%agemods])>0){
  fit.mv.bg<-eval(parse(
    text = paste("lm(ASO ~",
                 paste(uvmods.dc1m.bg.2,collapse = " + "),
                 "+Age+",
                 paste(
                   paste(uvmods.dc1m.bg.2[uvmods.dc1m.bg.2%in%agemods],"Age",sep = "*"),
                   collapse = " + "),
                 ", data = bg.dat)")
  ))
}else{
  fit.mv.bg<-eval(parse(
    text = paste("lm(ASO ~",
                 paste(uvmods.str.bg2,collapse = " + "),
                 ", data = bg.dat)")
  ))
}

fit.mv.bg.step<-stepAIC(fit.mv.bg,direction = "backward",trace = FALSE)
mvfit.dc.bg<-fit.mv.bg.step

restab<-as.data.frame(summary(fit.mv.bg.step)$coefficients)

fit.mv.bg.1m<-eval(parse(
  text = paste("lm(", 
               paste("ASO ~",
                     paste(uvmods.dc1m.bg.2,collapse = " + ")),
               ", data = bg.dat[bg.dat$Time==1,])")))

fit.mv.bg.1m.step<-stepAIC(fit.mv.bg.1m,direction = "backward",trace = FALSE)
mvfit.dc1m.bg<-fit.mv.bg.1m.step


restab<-as.data.frame(summary(fit.mv.bg.1m.step)$coefficients)
mvmods.dc1m.bg<-rownames(restab[restab$`Pr(>|t|)`<0.05,])
mvmods.dc1m.bg<-mvmods.dc1m.bg[!mvmods.dc1m.bg%in%"(Intercept)"]

fit.mv.bg.3m<-eval(parse(
  text = paste("lm(",
               paste("ASO ~",
                     paste(uvmods.dc3m.bg,collapse = " + ")),
               ", data = dc.dat[dc.dat$Time==3,])")))
mvfit.dc3m.bg<-fit.mv.bg.3m
sum.mv.bg.3m<-summary(fit.mv.bg.3m)
restab<-as.data.frame(sum.mv.bg.3m$coefficients)
mvmods.dc3m.bg<-rownames(restab[restab$`Pr(>|t|)`<0.05,])
mvmods.dc3m.bg<-mvmods.dc3m.bg[!mvmods.dc3m.bg%in%"(Intercept)"]



supMatTabs$titles$mv.str.bg<-paste0(
  "Striatum (matched samples) predictors of Thy1-haSyn (F = ",
  format(glance(mvfit.str.bg)$statistic,digits = 2),
  ", p = ",
  format(glance(mvfit.str.bg)$p.value,digits = 2),
  ")"
)
smdf<-summary(as.data.frame(report(mvfit.str.bg)))
smdf$Parameter<-gsub("str\\.","",smdf$Parameter)
supMatTabs$dfs$mv.str.bg<-smdf



supMatTabs$titles$mv.dc1m.bg<-paste0(
  "Colon (one month, matched samples) predictors of Thy1-haSyn (F = ",
  format(glance(mvfit.dc1m.bg)$statistic,digits = 2),
  ", p = ",
  format(glance(mvfit.dc1m.bg)$p.value,digits = 2),
  ")"
)
smdf<-summary(as.data.frame(report(mvfit.dc1m.bg)))
smdf$Parameter<-gsub("dc\\.","",smdf$Parameter)
supMatTabs$dfs$mv.dc1m.bg<-smdf


supMatTabs$titles$mv.dc3m.bg<-paste0(
  "Colon (three months, matched samples) predictors of Thy1-haSyn (F = ",
  format(glance(mvfit.dc3m.bg)$statistic,digits = 2),
  ", p = ",
  format(glance(mvfit.dc3m.bg)$p.value,digits = 2),
  ")"
)

smdf<-summary(as.data.frame(report(mvfit.dc3m.bg)))
smdf$Parameter<-gsub("dc\\.","",smdf$Parameter)
supMatTabs$dfs$mv.dc3m.bg<-smdf



# combined dc and str---------
univ.fit<-NULL

for(mod in names(str.MEs)){
  if(mod %in% agemods.str.bg){
    s.aso<-summary(lm(bg.dat$ASO~bg.dat[,mod]+bg.dat$Age+bg.dat[,mod]*bg.dat$Age))
  }else{
    s.aso<-summary(lm(bg.dat$ASO~bg.dat[,mod]))
    
  }
  s.aso1m<-summary(lm(bg.dat$ASO[bg.dat$Time==1]~bg.dat[bg.dat$Time==1,mod]))
  s.aso3m<-summary(lm(bg.dat$ASO[bg.dat$Time==3]~bg.dat[bg.dat$Time==3,mod]))
  
  univ.fit<-rbind(univ.fit,
                  c(
                    s.aso$coefficients[2,c("Estimate","Pr(>|t|)")],
                    s.aso1m$coefficients[2,c("Estimate","Pr(>|t|)")],
                    s.aso3m$coefficients[2,c("Estimate","Pr(>|t|)")]
                  ))
  
}
univ.fit<-as.data.frame(univ.fit)
names(univ.fit)<-c("asoEst","asoP","aso1mEst","aso1mP","aso3mEst","aso3mP")
rownames(univ.fit)<-names(str.MEs)
uvmods.str.bg<-rownames(univ.fit[univ.fit$asoP<0.05,])

uvmods.str1m.bg<-rownames(univ.fit[univ.fit$aso1mP<0.05,])
uvmods.str3m.bg<-rownames(univ.fit[univ.fit$aso3mP<0.05,])

bgmods.dc<-union(uvmods.dc.bg,union(uvmods.dc1m.bg,uvmods.dc3m.bg))
bgmods.dc<-bgmods.dc[!bgmods.dc%in%c("dc.cyan")]
bgmods.str<-union(uvmods.str.bg,union(uvmods.str1m.bg,uvmods.str3m.bg))

agemods.bg<-union(agemods.dc.bg,agemods.str.bg)
agemods.bg<-union(agemods.dc,agemods.str)

bgmods<-union(uvmods.dc.bg,uvmods.str.bg)

if(length(bgmods[bgmods%in%agemods.bg])>0){
  bgfit.comb<-eval(parse(
    text = paste("lm(ASO ~",
                 paste(bgmods,collapse = " + "),
                 "+Age+",
                 paste(
                   paste(bgmods[bgmods%in%agemods.bg],"Age",sep = "*"),
                   collapse = " + "),
                 ", data = bg.dat)")
  ))
}else{
  bgfit.comb<-eval(parse(
    text = paste("lm(ASO ~",
                 paste(bgmods,collapse = " + "),
                 ", data = bg.dat)")
  ))
}
summary(bgfit.comb)

bgfit.comb<-stepAIC(bgfit.comb,direction = "backward",trace = FALSE)
summary(bgfit.comb)



bgfit.dc<-lm(ASO~dc.yellow + Age + dc.yellow:Age, data = bg.dat)
summary(bgfit.dc)

bgfit.str<-lm(ASO~str.royalblue + str.red + str.cyan +Age, data = bg.dat)

summary(bgfit.str)

aov.bgfit<-anova(bgfit.str,bgfit.comb)

bgfit1m.comb<-eval(parse(
  text = paste("lm(", 
               paste("ASO ~",
                     paste(mvmods.str.bg,collapse = " + "),
                     " + ",
                     paste(mvmods.dc1m.bg,collapse = " + ")),
               ", data = bg.dat[bg.dat$Time==1,])")))
bgfit1m.comb<-stepAIC(bgfit1m.comb,direction = "backward",trace = FALSE)
summary(bgfit1m.comb)
bgfit.dc1m<-lm(ASO ~dc.darkmagenta, data = bg.dat[bg.dat$Time==1,])
bgfit.str1m<-lm(ASO ~str.red, data = bg.dat[bg.dat$Time==1,])

aov.bgfit1m<-anova(bgfit.str1m,bgfit1m.comb)







save(bgfit.comb,bgfit.dc,bgfit.str,
     bgfit.dc1m,bgfit.str1m,bgfit1m.comb,
     mvfit.dc, mvfit.dc1m, mvfit.dc3m,
     mvfit.str, 
     supMatTabs,
     file = "regResults.rda")

saveRDS(supMatTabs,file = "SupMatRegTabs.rds")

save(agemods.dc, agemods.dc.bg,
     agemods.str, agemods.str.bg,
     bgmods.dc,bgmods.str,
     mvmods.dc1m,mvmods.dc3m,
     uvmods.dc, uvmods.dc.bg,
     uvmods.dc1m, uvmods.dc1m.bg,
     uvmods.dc3m, uvmods.dc3m.bg,
     uvmods.str, uvmods.str.bg,
     file = "regressionModlists.rda")


