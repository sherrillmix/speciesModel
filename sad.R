library(parallel)
#radfit likes to fail on madelbrot fit
source('rad.zipfbrot.R')


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
table(apply(aics[isEnough,],1,which.min))
table(apply(aics2[isEnough2,],1,which.min))
mandels<-do.call(rbind,lapply(rads,function(xx)xx$model[[5]]$coef))
mandels2<-do.call(rbind,lapply(rads2,function(xx)xx$model[[5]]$coef))
zipf<-do.call(rbind,lapply(rads,function(xx)xx$model[[4]]$coef))
zipf2<-do.call(rbind,lapply(rads2,function(xx)xx$model[[4]]$coef))
rads[[50]]
tapply(mandels[isEnough,2],info[rownames(mandels)[isEnough],'class'],mean)


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
dev.off()

