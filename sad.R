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


rads<-mclapply(speciesAbund,radfit,mc.cores=20)
rads2<-mclapply(speciesAbund2,radfit,mc.cores=20)
rads[isEnough][1]
rads2[isEnough2][1]
aics<-do.call(rbind,lapply(rads,function(xx)structure(sapply(xx$models,function(yy)yy$aic),.Names=sapply(xx$models,function(yy)yy$model))))
aics2<-do.call(rbind,lapply(rads2,function(xx)structure(sapply(xx$models,function(yy)yy$aic),.Names=sapply(xx$models,function(yy)yy$model))))
table(apply(aics[isEnough,],1,which.min),info[rownames(aics[isEnough,]),'class'])
table(apply(aics2[isEnough2,],1,which.min))
mandels<-do.call(rbind,lapply(rads,function(xx)xx$model[[5]]$coef))
mandels2<-do.call(rbind,lapply(rads2,function(xx)xx$model[[5]]$coef))
zipf<-do.call(rbind,lapply(rads,function(xx)xx$model[[4]]$coef))
zipf2<-do.call(rbind,lapply(rads2,function(xx)xx$model[[4]]$coef))
rads[[50]]
tapply(mandels[isEnough,2],info[rownames(mandels)[isEnough],'class'],mean)
pbs<-lapply(speciesAbund,fitpowbend)

pows<-mclapply(speciesAbund[isEnough],fitpower,mc.cores=20)
pows2<-mclapply(speciesAbund2[isEnough2],fitpower,mc.cores=20)
powbends<-mclapply(speciesAbund[isEnough],fitpowbend,mc.cores=20)
powbends2<-mclapply(speciesAbund[isEnough2],fitpowbend,mc.cores=20)
pdf('out/neutralTest.pdf')
plot(pows[[1]])
plot(fitpowbend(speciesAbund[isEnough][[1]]))
plot(fitmzsm(speciesAbund[isEnough][[1]]))
plot(fitvolkov(speciesAbund[isEnough][[1]]))
dev.off()

powCoef<-sapply(pows,function(xx)xx@coef)
powCoef2<-sapply(pows2,function(xx)xx@coef)
names(powCoef2)<-sub('.s$','',names(powCoef2))
powCoef2<-powCoef2[info2[names(pows2),'study']!='bushman'&names(pows2) %in% rownames(info2)]

fitters<-c('Broken stick'=fitbs2,'Log series'=fitls2,'Neutral (MZSM)'=fitmzsm2,'Power'=fitpower,'Power bend'=fitpowbend2,'Poisson lognormal'=fitpoilog,'Negative binom'=fitnbinom2,'Geometric'=fitgeom2,'Log normal'=fitlnorm,'Gamma'=fitgamma2,'Weibull'=fitweibull)
#fitters<-c('Zipf'=fitzipf,'Zipf-Mandelbrot'=fitmand,'Rank broken stick'=fitrbs,'Geometric'=fitgs) #,fitvolkov
fits<-mclapply(speciesAbund[isEnough],function(xx){
  cat('.')
  cacheFile<-sprintf('cache/%s.Rdat',digest(xx))
  if(TRUE || !exists(cacheFile)){
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
    else return(BIC(yy))
  })
}))
bicDiff<-t(apply(bics,1,function(xx)xx-min(xx,na.rm=TRUE)))
bicDiff[is.infinite(bicDiff)]<-NA
bicDiff<-bicDiff[do.call(order,info[rownames(bics),c('phylum','supersuperclass','superclass','class','superorder','clade','order','family','genus','species')]),]
#bicDiff[is.na(bicDiff)]<-max(bicDiff+1,na.rm=TRUE)
colnames(bicDiff)<-names(fitters)
bicCondense<-do.call(rbind,by(bicDiff,info[rownames(bicDiff),'common'],function(xx)apply(xx,2,mean,na.rm=TRUE)))
bicCondense<-bicCondense[orderIn(rownames(bicCondense),info[rownames(bicDiff),'common']),]
pdf('out/picks.pdf',height=10)
  par(mar=c(7,10,1,1))
  image(1:ncol(bicCondense),1:nrow(bicCondense),t(log2(1+bicCondense)),xaxt='n',col=colorRampPalette(c('red','blue'))(100),ylab='',xlab='',yaxt='n')
  axis(2,1:nrow(bicCondense),rownames(bicCondense),las=2,cex.axis=.7)
  axis(1,1:ncol(bicCondense),colnames(bicCondense),las=2)
dev.off()

#color brewer
fitCols<-c('#8dd3c7CC','#ffffb3CC','#bebadaCC','#fb8072CC','#80b1d3CC','#fdb462CC','#b3de69CC','#fccde5CC','#bc80bdCC')
names(fitCols)<-c('Broken stick','Log norm','Log series','Neutral (MZSM)','Neutral (Volkov)','Power series','Power bend','Gamma','Neg binom')
targets<-c('Stinkbug','Human','Goat','Right Whale')
ids<-sapply(targets,function(xx)which(info[rownames(bicDiff),'common']==xx)[1])
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

preds<-mclapply(targets,function(target){
  ii<-ids[target]
  lapply(fits[[ii]],radpred)
},mc.cores=5)
names(preds)<-targets

pdf('out/fitPicks.pdf')
  layout(matrix(c(rep(5,4),1:4),ncol=2),width=c(.7,.3))
  for(target in targets){
    ii<-ids[target]
    thisRad<-rad(speciesAbund[isEnough][[ii]])
    plot(thisRad$rank,thisRad$abund,main=target,las=1,log='y',xlab='OTU rank',ylab='OTU abundance')
    isNA<-sapply(1:length(fits[[ii]]),function(xx){
      if(is.null(fits[[ii]][[xx]]))return(TRUE)
      lines(preds[[target]][[xx]],col=fitCols[xx],lwd=2)
      return(FALSE)
    })
    legend('topright',names(fitCols)[!isNA],col=fitCols[!isNA],lwd=2)
  }
  par(mar=c(7,10,1,1))
  image(1:ncol(bicCondense),1:nrow(bicCondense),t(log2(1+bicCondense)),xaxt='n',col=colorRampPalette(c('red','blue'))(100),ylab='',xlab='',yaxt='n')
  axis(2,1:nrow(bicCondense),rownames(bicCondense),las=2,cex.axis=.7)
  axis(1,1:ncol(bicCondense),colnames(bicCondense),las=2)
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


#http://andrewgelman.com/2016/06/11/log-sum-of-exponentials/
log_sum_exp<-function(u, v) max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))

pZM<-function(N,q,s){
  #fs<-1/(1:N+q)^s
  logFs<- -s*log(1:N+q)
  H<-Reduce(log_sum_exp,logFs)
  exp(logFs-H)
}
pRadZM<-function(ns,q,s,log=TRUE){
  if(q<=-1)return(ifelse(log,-Inf,0))
  N<-sum(ns)
  rad<-table(factor(ns,levels=1:N))
  ps<-pZM(N,q,s)
  out<-dmultinom(rad,sum(rad),ps,log=log)
  return(out)
}
fitZM<-function(ns){
  fit<-withCallingHandlers(
    nlm(function(qs,ns)-pRadZM(ns,qs[1],qs[2]),c(0,1),ns),
    warning=function(w)if(grepl('Inf replaced by maximum positive value',w))invokeRestart('muffleWarning')
  )
  return(c('q'=fit$estimate[1],'s'=fit$estimate[2],'log-likelihood'=-fit$minimum))
}

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
