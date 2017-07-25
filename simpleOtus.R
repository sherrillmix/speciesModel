if(!exists('combo'))source('readCombo.R')

pdf('out/simple.pdf')
  for(ii in c('all',unique(combo$study))){
    if(ii=='all')thisCombo<-combo
    else thisCombo<-combo[combo$study==ii,]
    plot(thisCombo$weight,thisCombo$otu,log='xy',xaxt='n',yaxt='n',xlab='Weight (kg)',ylab='Rarefied OTUs',main=ii,pch=21,bg='#00000099',col=NA)
    logAxis(1,addExtra=TRUE)
    logAxis(2,las=1)
    y<-log10(thisCombo$otu)
    x<-log10(thisCombo$weight)
    tmp<-data.frame(x,y)
    fit<-lm(y~x,dat=tmp)
    fakeWeight<-seq(min(thisCombo$log.weight),max(thisCombo$log.weight),length.out=100)
    pred<-predict(fit,data.frame(x=fakeWeight),interval='confidence')
    lines(10^fakeWeight,10^pred[,'fit'],col='#FF000099')
    polygon(c(10^fakeWeight,rev(10^fakeWeight)),c(10^pred[,'lwr'],rev(10^pred[,'upr'])),col='#FF000033',border=NA)
  }
dev.off()
