library(parallel)
library(sads)
library(vegan)
library(digest)
library(dnar)
#radfit likes to fail on madelbrot fit
source('rad.zipfbrot.R')
source('functions.R')


if(!file.exists('work/rareN.csv'))source('readData.R')
load('work/speciesAbund.Rdat')
load('work/speciesAbund2.Rdat')
info<-read.csv('data/islandGut discovery metadata - map.tsv.csv',row.names=1,stringsAsFactors=FALSE)
info2<-read.csv('data/islandGut validation metadata - map.tsv.csv',row.names=1,stringsAsFactors=FALSE)

nReads<-100
n<-sapply(speciesAbund,sum)
isEnough<-n>nReads
n2<-sapply(speciesAbund2,sum)
isEnough2<-n2>nReads


fitters<-c('Broken stick'=fitbs2,'Geometric'=fitgeom2,'Log series'=fitls2,'Neutral (MZSM)'=fitmzsm2,'Power'=fitpower,'Power bend'=fitpowbend2,'Log normal'=fitlnorm,'Poisson lognormal'=fitpoilog2,'Gamma'=fitgamma2,'Weibull'=fitweibull) #,'Negative binom'=fitnbinom2
load('work/deblurAbund.Rdat')
#fitters<-c('Zipf'=fitzipf,'Zipf-Mandelbrot'=fitmand,'Rank broken stick'=fitrbs,'Geometric'=fitgs) #,fitvolkov
fits<-mclapply(speciesAbund[isEnough],function(xx){
  cat('.')
  cacheFile<-sprintf('cache/%s.Rdat',digest(xx))
  if(!exists(cacheFile)){
    out<-lapply(fitters,function(func,xx){
      tryCatch(func(xx),error=function(e)return(NULL))
    },xx)
    save(out,file=cacheFile)
  }else{
    tmp<-new.env()
    load(cacheFile,tmp)
    out<-get('out',tmp)
  }
  return(out)
},mc.cores=30,mc.preschedule=FALSE)

bics<-do.call(rbind,lapply(fits,function(xx){
  sapply(xx,function(yy){
    if(is.null(yy))return(NA)
    else return(AIC(yy))
  })
}))
bicDiff<-t(apply(bics,1,function(xx)xx-min(xx,na.rm=TRUE)))
bicDiff[is.infinite(bicDiff)]<-NA
bicDiff<-bicDiff[do.call(order,cbind(info[rownames(bics),c('phylum','supersuperclass','superclass','class','superorder','clade','order','family','genus','species')])),]
#bicDiff[is.na(bicDiff)]<-max(bicDiff+1,na.rm=TRUE)
colnames(bicDiff)<-names(fitters)
bicCondense<-do.call(rbind,by(bicDiff,info[rownames(bicDiff),'common'],function(xx)apply(xx,2,mean,na.rm=TRUE)))
bicCondense<-bicCondense[orderIn(rownames(bicCondense),info[rownames(bicDiff),'common']),]
rownames(bicCondense)<-sub('TigerSha','Tiger sha',sub('SSBShark','Sandbar shark',sub('SBUShark','Bull shark',sub('Vietnames ','Vietnamese ',rownames(bicCondense)))))
bicCondense<-bicCondense[rownames(bicCondense)!='unknown_or_none',]
pdf('out/picks.pdf',height=10)
  par(mar=c(7,10,1,1))
  image(1:ncol(bicCondense),1:nrow(bicCondense),t(log2(1+bicCondense)),xaxt='n',col=colorRampPalette(c('red','blue'))(100),ylab='',xlab='',yaxt='n')
  axis(2,1:nrow(bicCondense),rownames(bicCondense),las=2,cex.axis=.7)
  axis(1,1:ncol(bicCondense),colnames(bicCondense),las=2)
dev.off()

#color brewer
#fitCols<-c('#8dd3c7CC','#ffffb3CC','#bebadaCC','#fb8072CC','#80b1d3CC','#fdb462CC','#b3de69CC','#fccde5CC','#bc80bdCC')
fitCols<-c('#a6cee3AA','#1f78b4AA','#b2df8aAA','#33a02cAA','#fb9a99AA','#e31a1cAA','#fdbf6fAA','#ff7f00AA','#cab2d6AA','#6a3d9aAA','#ffff99AA')[1:ncol(bicDiff)]
names(fitCols)<-colnames(bicDiff)
#c('Broken stick','Log norm','Log series','Neutral (MZSM)','Neutral (Volkov)','Power series','Power bend','Gamma','Neg binom')
targets<-c('Stinkbug'='Stinkbug','Human'='Human','TigerShark'='Tiger shark','Right Whale'='Right Whale')
ids<-sapply(names(targets),function(xx)which(info[rownames(bicDiff),'common']==xx)[1])
names(ids)<-targets
pdf('out/fits.pdf')
for(target in targets){
  ii<-ids[target]
  thisRad<-rad(speciesAbund[isEnough][[ii]])
  plot(thisRad$rank,thisRad$abund,main=target,las=1,log='y',xlab='OTU rank',ylab='OTU abundance')
  isNA<-sapply(1:length(fits[[ii]]),function(xx){
    if(is.null(fits[[ii]][[xx]]))return(TRUE)
    lines(radpred(fits[[ii]][[xx]]),col=fitCols[xx],lwd=2)
    return(FALSE)
  })
  legend('topright',names(fitCols)[!isNA],col=fitCols[!isNA],lwd=2)
}
dev.off()

preds<-mclapply(targets,function(target,...){
  ii<-ids[target]
  message(target)
  out<-mclapply(names(fits[[ii]]),function(jj){
    message(jj)
    lazyRadPred(fits[[ii]][[jj]])
  },mc.cores=10)
  names(out)<-names(fits[[ii]])
  return(out)
},mc.cores=5)
names(preds)<-targets

pdf('out/fitPicks.pdf',width=8,height=10)
  layout(matrix(c(rep(5,4),1:4),ncol=2),width=c(.7,.3))
  par(mar=c(3,3.5,1,1))
  targetTops<-targetBottoms<-targetLefts<-c()
  legendAdded<-FALSE
  for(target in targets[orderIn(targets,rownames(bicCondense),decreasing=TRUE)]){
    ii<-ids[target]
    print(names(speciesAbund[isEnough])[[ii]])
    thisRad<-rad(speciesAbund[isEnough][[ii]])
    plot(thisRad$rank,thisRad$abund,main=target,las=1,log='y',xlab='',ylab='OTU abundance',mgp=c(2.5,.7,0))
    title(xlab='OTU rank',mgp=c(1.6,1,0))
    isNA<-sapply(names(fitCols),function(xx){
      if(is.null(fits[[ii]][[xx]]))return(TRUE)
      thisPred<-preds[[target]][[xx]]
      thisPred[is.infinite(thisPred[,'abund']),'abund']<-sum(thisRad)*2
      lines(thisPred,col=fitCols[xx],lwd=2)
      return(FALSE)
    })
    targetLefts[target]<-grconvertX(par('usr')[1],to='ndc')
    targetTops[target]<-grconvertY(10^par('usr')[4],to='ndc')
    targetBottoms[target]<-grconvertY(10^par('usr')[3],to='ndc')
    if(!legendAdded)legend('topright',names(fitCols)[!isNA],col=fitCols[!isNA],lwd=2,bty='n',cex=.95)
    legendAdded<-TRUE
  }
  par(mar=c(7,10,1,1.5))
  cols<-colorRampPalette(c('red','blue'))(200)
  breaks<-seq(min(log10(1+bicCondense),na.rm=TRUE),max(log10(1+bicCondense),na.rm=TRUE),length.out=201)
  image(1:ncol(bicCondense),1:nrow(bicCondense),t(log10(1+bicCondense)),xaxt='n',col=cols,breaks=breaks,ylab='',xlab='',yaxt='n')
  box()
  axis(2,1:nrow(bicCondense),rownames(bicCondense),las=2,cex.axis=.7)
  slantAxis(1,1:ncol(bicCondense),colnames(bicCondense),srt=-30,adj=c(.1,.5))
  for(target in targets){
    y1<-which(rownames(bicCondense)==target)+c(.5,-.5)
    x1<-rep(par('usr')[2],2)
    x2<-rep(grconvertX(targetLefts[target],'ndc','user'),2)
    y2<-grconvertY(c(targetTops[target],targetBottoms[target]),'ndc','user')
    #polygon(c(x1,rev(x2)),c(y1,rev(y2)),border='#00000099',col=NA,xpd=NA,lty=2)
    polygon(c(x1,rev(x2)),c(y1,rev(y2)),border='#00000011',col='#00000011',xpd=NA)
  }
  insetScale(breaks,cols,c(.015,.015,.025,.3),main='Difference from minimum AIC',at=log10(c(0,10^(1:5))+1),labels=c(0,sapply(1:5,function(xx)as.expression(bquote(10^.(xx))))),offset=.01)
dev.off()

pdf('out/params.pdf')
vpPlot(info[rownames(zipf[isEnough,]),'class'],log10(-zipf[isEnough,'gamma']),las=2,ylab='gamma parameter Zipf')
vpPlot(info[rownames(mandels[isEnough&!is.infinite(mandels[,'c']),]),'class'],log10(-mandels[isEnough&!is.infinite(mandels[,'c']),'gamma']),las=2,ylab='gamma parameter Zipf-Mandelbrot')
vpPlot(info[rownames(mandels[isEnough&!is.infinite(mandels[,'c']),]),'class'],mandels[isEnough&!is.infinite(mandels[,'c']),'beta'],las=2,ylab='beta parameter Zipf-Mandelbrot')
ylim<-range(c(powCoef,powCoef2))
vpPlot(info[names(pows),'class'],powCoef,las=2,ylab='Zipf s',ylim=ylim)
studyCols<-rainbow(length(unique(info2[names(powCoef2),'study'])))
names(studyCols)<-unique(info2[names(powCoef2),'study'])
vpPlot(info2[names(powCoef2),'class'],powCoef2,las=2,ylab='Zipf s',ylim=ylim,col=studyCols[info2[names(powCoef2),'study']])
legend('topleft',names(studyCols),col=studyCols,pch=1)
dev.off()


zms<-do.call(rbind,mclapply(speciesAbund[isEnough],fitZM,mc.cores=20))

library(vipor)
pdf('test.pdf')
ylim<-range(log10(zms[,'s']))
ylim[1]<-min(ylim[1],0)
vpPlot(info[rownames(zms),'class'],log10(zms[,'s']),las=2,ylab='s parameter Zipf-Mandelbrot',yaxt='n',ylim=ylim)
prettyY<-c(1,5,10,20)
axis(2,log10(prettyY),prettyY,,las=1)
vpPlot(info[rownames(zms),'class'],zms[,'q'],las=2,ylab='q parameter Zipf-Mandelbrot')
dev.off()

pdf('fits.pdf')
#for(ii in 1:sum(isEnough)){
for(ii in 1:10){
  counts<-as.vector(table(factor(speciesAbund[isEnough][[ii]],levels=1:max(speciesAbund[isEnough][[ii]]))))
  preds<-pZM(sum(speciesAbund[isEnough][[ii]]),zms[ii,'q'],zms[ii,'s'])
  plot((counts+1)/sum(counts+1),log='xy',ylim=range(preds),main=names(isEnough)[ii])
  points(preds,pch='-',col='red')
}
dev.off()

#deblur
load('work/deblurAbund.Rdat')

nDeblur<-sapply(deblurAbund,sum)
isEnoughDeblur<-nDeblur>nReads
nDeblur2<-sapply(deblurAbund2,sum)
isEnoughDeblur2<-nDeblur2>nReads
nDeblur4<-sapply(deblurAbund4,sum)
isEnoughDeblur4<-nDeblur4>nReads
deblurFits<-mclapply(deblurAbund[isEnoughDeblur],function(xx){
  cat('.')
  out<-lapply(fitters,function(func,xx){
    tryCatch(func(xx),error=function(e)return(NULL))
  },xx)
  return(out)
},mc.cores=30,mc.preschedule=FALSE)
deblurFits2<-mclapply(deblurAbund2[isEnoughDeblur2],function(xx){
  cat('.')
  out<-lapply(fitters,function(func,xx){
    tryCatch(func(xx),error=function(e)return(NULL))
  },xx)
  return(out)
},mc.cores=30,mc.preschedule=FALSE)
deblurFits4<-mclapply(deblurAbund4[isEnoughDeblur4],function(xx){
  cat('.')
  out<-lapply(fitters,function(func,xx){
    tryCatch(func(xx),error=function(e)return(NULL))
  },xx)
  return(out)
},mc.cores=30,mc.preschedule=FALSE)

plotBics<-function(fits,speciesAbund,info,outFile){
  bics<-do.call(rbind,lapply(fits,function(xx){
    sapply(xx,function(yy){
      if(is.null(yy))return(NA)
      else return(AIC(yy))
    })
  }))
  bicDiff<-t(apply(bics,1,function(xx)xx-min(xx,na.rm=TRUE)))
  bicDiff[is.infinite(bicDiff)]<-NA
  bicDiff<-bicDiff[do.call(order,cbind(info[rownames(bics),c('phylum','supersuperclass','superclass','class','superorder','clade','order','family','genus','species')])),]
  colnames(bicDiff)<-names(fitters)
  #need to reorder to match bicDiff reordering
  speciesAbund<-speciesAbund[rownames(bicDiff)]
  fits<-fits[rownames(bicDiff)]
  bicCondense<-do.call(rbind,by(bicDiff,info[rownames(bicDiff),'common'],function(xx)apply(xx,2,mean,na.rm=TRUE)))
  bicCondense<-bicCondense[orderIn(rownames(bicCondense),info[rownames(bicDiff),'common']),]
  rownames(bicCondense)<-sub('TigerSha','Tiger sha',sub('SSBShark','Sandbar shark',sub('SBUShark','Bull shark',sub('Vietnames ','Vietnamese ',rownames(bicCondense)))))
  bicCondense<-bicCondense[rownames(bicCondense)!='unknown_or_none',]
  #color brewer
  fitCols<-c('#a6cee3AA','#1f78b4AA','#b2df8aAA','#33a02cAA','#fb9a99AA','#e31a1cAA','#fdbf6fAA','#ff7f00AA','#cab2d6AA','#6a3d9aAA','#ffff99AA')[1:ncol(bicDiff)]
  names(fitCols)<-colnames(bicDiff)
  targets<-c('Stinkbug'='Stinkbug','Human'='Human','TigerShark'='Tiger shark','Right Whale'='Right Whale')
  ids<-sapply(names(targets),function(xx)which(info[rownames(bicDiff),'common']==xx)[1])
  names(ids)<-targets
  preds<-mclapply(targets,function(target,...){
    ii<-ids[target]
    message(target)
    out<-mclapply(names(fits[[ii]]),function(jj){
      message(jj)
      #time out after 1 hour
      tryCatch(R.utils::withTimeout(lazyRadPred(fits[[ii]][[jj]]),timeout=3600),TimeoutException=function(ex){warning('Time out ',jj);NULL})
    },mc.cores=10)
    names(out)<-names(fits[[ii]])
    return(out)
  },mc.cores=5)
  names(preds)<-targets
  pdf(outFile,width=8,height=10)
    layout(matrix(c(rep(5,4),1:4),ncol=2),width=c(.7,.3))
    par(mar=c(3,3.5,1,1))
    targetTops<-targetBottoms<-targetLefts<-c()
    legendAdded<-FALSE
    for(target in targets[orderIn(targets,rownames(bicCondense),decreasing=TRUE)]){
      ii<-ids[target]
      thisRad<-rad(speciesAbund[[ii]])
      plot(thisRad$rank,thisRad$abund,main=target,las=1,log='y',xlab='',ylab='OTU abundance',mgp=c(2.5,.7,0),yaxt='n')
      logAxis(2,las=1)
      title(xlab='OTU rank',mgp=c(1.6,1,0))
      isNA<-sapply(names(fitCols),function(xx){
        if(is.null(fits[[ii]][[xx]]))return(TRUE)
        if(is.null(preds[[target]][[xx]]))return(TRUE)
        thisPred<-preds[[target]][[xx]]
        thisPred[is.infinite(thisPred[,'abund']),'abund']<-sum(thisRad)*2
        lines(thisPred,col=fitCols[xx],lwd=2)
        return(FALSE)
      })
      targetLefts[target]<-grconvertX(par('usr')[1],to='ndc')
      targetTops[target]<-grconvertY(10^par('usr')[4],to='ndc')
      targetBottoms[target]<-grconvertY(10^par('usr')[3],to='ndc')
      if(!legendAdded)legend('topright',names(fitCols)[!isNA],col=fitCols[!isNA],lwd=2,bty='n',cex=.95)
      legendAdded<-TRUE
    }
    par(mar=c(7,10,1,1.5))
    cols<-colorRampPalette(c('red','blue'))(200)
    breaks<-seq(min(log10(1+bicCondense),na.rm=TRUE),max(log10(1+bicCondense),na.rm=TRUE),length.out=201)
    image(1:ncol(bicCondense),1:nrow(bicCondense),t(log10(1+bicCondense)),xaxt='n',col=cols,breaks=breaks,ylab='',xlab='',yaxt='n')
    box()
    axis(2,1:nrow(bicCondense),rownames(bicCondense),las=2,cex.axis=.7)
    slantAxis(1,1:ncol(bicCondense),colnames(bicCondense),srt=-30,adj=c(.1,.5))
    for(target in targets){
      y1<-which(rownames(bicCondense)==target)+c(.5,-.5)
      x1<-rep(par('usr')[2],2)
      x2<-rep(grconvertX(targetLefts[target],'ndc','user'),2)
      y2<-grconvertY(c(targetTops[target],targetBottoms[target]),'ndc','user')
      #polygon(c(x1,rev(x2)),c(y1,rev(y2)),border='#00000099',col=NA,xpd=NA,lty=2)
      polygon(c(x1,rev(x2)),c(y1,rev(y2)),border='#00000011',col='#00000011',xpd=NA)
    }
    ticks<-1:floor(max(log10(bicCondense),na.rm=TRUE))
    insetScale(breaks,cols,c(.015,.015,.025,.3),main='Difference from minimum AIC',at=log10(c(0,10^(ticks)+1)),labels=c(0,sapply(ticks,function(xx)as.expression(bquote(10^.(xx))))),offset=.01)
  dev.off()
}

plotBics(deblurFits,deblurAbund,info,'out/deblurFit.pdf')
plotBics(deblurFits2,deblurAbund2,info,'out/deblur2Fit.pdf')
plotBics(deblurFits4,deblurAbund4,info,'out/deblur4Fit.pdf')

#dada2
load('work/dadaAbund.Rdat')
nDada<-sapply(dadaAbund,sum)
isEnoughDada<-nDada>nReads
dadaFits<-mclapply(dadaAbund[isEnoughDada],function(xx){
  cat('.')
  xx<-xx[xx>1]
  out<-lapply(fitters,function(func,xx){
    tryCatch(func(xx,trunc=1),error=function(e)return(NULL))
  },xx)
  return(out)
},mc.cores=30,mc.preschedule=FALSE)

plotBics(dadaFits,lapply(dadaAbund,function(xx)xx[xx>1]),info,'out/dadaFit.pdf')



