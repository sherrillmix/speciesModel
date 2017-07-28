library(dnar)

if(!exists('combo'))source('readCombo.R')

lms<-list('sample'=list(),'species'=list())
pdf('out/simple.pdf')
  for(ii in c('all','bushman','others',unique(combo$study[combo$study!='bushman']))){
    if(ii=='all')thisCombo<-combo
    else if(ii=='others')thisCombo<-combo[combo$study!='bushman',]
    else thisCombo<-combo[combo$study==ii,]
    thisAves<-data.frame('otu'=tapply(thisCombo$otu,thisCombo$species,function(xx)exp(mean(log(xx)))),'weight'=tapply(thisCombo$weight,thisCombo$species,function(xx)exp(mean(log(xx)))))
    thisAves$log.weight<-log10(thisAves$weight)
    thisAves$log.otu<-log10(thisAves$otu)
    thisSets<-list('sample'=thisCombo,'species'=thisAves)
    for(jj in names(thisSets)){
      thisSet<-thisSets[[jj]]
      plot(thisSet$weight,thisSet$otu,log='xy',xaxt='n',yaxt='n',xlab='Weight (kg)',ylab='Rarefied OTUs',main=sprintf('%s %s',ii,jj),pch=21,bg='#00000099',col=NA,xlim=range(combo$weight),ylim=range(combo$otu))
      logAxis(1,addExtra=TRUE)
      logAxis(2,las=1)
      otu<-log10(thisSet$otu)
      weight<-log10(thisSet$weight)
      tmp<-data.frame(otu,weight)
      fit<-lm(otu~weight,dat=tmp)
      lms[[jj]][[ii]]<-fit
      fakeWeight<-seq(min(thisSet$log.weight),max(thisSet$log.weight),length.out=100)
      pred<-predict(fit,data.frame(weight=fakeWeight),interval='confidence')
      lines(10^fakeWeight,10^pred[,'fit'],col='#FF000099')
      polygon(c(10^fakeWeight,rev(10^fakeWeight)),c(10^pred[,'lwr'],rev(10^pred[,'upr'])),col='#FF000033',border=NA)
    }
  }
dev.off()

write.csv(do.call(rbind,lapply(lms[['sample']],function(xx)summary(xx)$coefficients['weight',])),'out/slopesSample.csv')
write.csv(do.call(rbind,lapply(lms[['species']],function(xx)summary(xx)$coefficients['weight',])),'out/slopesSpecies.csv')
