library(rstan)
library(vipor)
library(dnar)
library(vioplot)

#use multiple cores (careful with big models on small computers)
options(mc.cores = parallel::detectCores())
#cache compilations (doesn't actually work that well, at least with model_code=)
rstan_options(auto_write = TRUE)


stanCode<-'
  data{
    int<lower=0> n;
    real otus[n]; //log
    int<lower=1> ids[n];
    int<lower=0> nIds;
    vector[nIds] weights; //log
    int<lower=1> classes[nIds];
    int<lower=0> nClasses;
  }
  parameters {
    vector[nClasses] intercepts;
    real<lower=0> metaSigma;
    //vector<lower=0>[nClasses] sigmas;
    real<lower=0> sigma;
    vector[nIds] rawMus;
    vector[nClasses] rawBetas;
    real metaBeta;
    real<lower=0> metaBetaSd;
    //real<lower=0> speciesVar;
  }
  transformed parameters {
    vector[nIds] mus;
    vector[nClasses] betas;
    betas=metaBeta+metaBetaSd*rawBetas;
    mus = intercepts[classes] + betas[classes] .* weights + rawMus*metaSigma;
  }
  model {
    //speciesVar ~ gamma(1,3);
    //sigmas ~ gamma(1,speciesVar);
    //sigmas ~ gamma(1,3);
    sigma ~ gamma(1,.1);
    rawBetas ~ normal(0,1);
    metaBetaSd ~ gamma(1,.1);
    metaSigma ~ gamma(1,3);
    rawMus ~ normal(0,1);
    //otus ~ normal(mus[ids],sigmas[classes[ids]]);
    otus ~ normal(mus[ids],sigma);
  }
'



dat<-read.csv('tableScott.discovery.csv',row.names=1,stringsAsFactors=FALSE)
colnames(dat)[colnames(dat)=='host_species_common_name']<-'common'
dat$study<-'bushman'
#CAREFUL excluding single sample classes for now
#dat<-dat[ave(dat$class,dat$class,FUN=length)>1,]
speciesIds<-1:length(unique(dat$common))
names(speciesIds)<-unique(dat$common)
dat$speciesId<-speciesIds[dat$common]
dat$group<-dat$class
classIds<-1:length(unique(dat$class))
names(classIds)<-unique(dat$class)
dat$classId<-classIds[dat$class]
speciesIdClasses<-sapply(names(speciesIds),function(xx)dat[dat$common==xx,][1,'class'])
dat$taxa<-tolower(apply(dat[,c('phylum','class','order','family','genus','species')],1,paste,collapse=' '))


plot(dat$log.weight,dat$log.otu,pch=21,bg=rainbow(max(dat$speciesId))[dat$speciesId],col=rainbow(max(dat$classId))[dat$classId],lwd=2)
meanWeight<-tapply(dat$log.weight,dat$common,mean)
meanOtu<-tapply(dat$log.otu,dat$common,mean)
plot(meanWeight, meanOtu)
summary(lm(meanOtu~meanWeight))

fit<-stan(
  model_code=stanCode,
  data=list(
    otus=dat$log.otu,
    n=nrow(dat),
    ids=dat$speciesId,
    nIds=max(dat$speciesId),
    weights=tapply(dat$log.weight,dat$speciesId,mean),
    classes=tapply(dat$classId,dat$speciesId,'[[',1),
    nClasses=max(dat$classId)
  ),
  control=list(adapt_delta=.99),
  chains=32,iter=100000,thin=100
)
print(fit)

plot(fit)
print(traceplot(fit))
pdf('trace.pdf');print(traceplot(fit,par=c('metaBeta','betas','metaSigma','sigma')));plot(fit,par=c('metaBeta','betas'));dev.off()
sims<-extract(fit)
hist(sims[['metaBeta']])
speciesMeans<-apply(sims[['mus']],2,mean)
weights<-tapply(dat$log.weight,dat$speciesId,mean)
speciesCols<-rainbow(max(dat$speciesId))
xlim<-range(dat$log.weight)
ylim<-range(dat$log.otu)
pdf('data.pdf',width=12,height=12)
plot(dat$log.weight,dat$log.otu,pch=21,bg=speciesCols[dat$speciesId],col=rainbow(max(dat$classId))[dat$classId],lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',xlim=xlim,ylim=ylim)
segments(weights-.2,speciesMeans,weights+.2,speciesMeans,col=speciesCols,lwd=2)
for(ii in 1:max(dat$classId)){
  withAs(dat=dat[dat$classId==ii,],plot(dat$log.weight,dat$log.otu,pch=21,bg=speciesCols[dat$speciesId],col=NA,lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',main=names(classIds)[ii],xlim=xlim,ylim=ylim,cex=2))
  classSelect<-speciesIdClasses==names(classIds)[ii]
  segments(weights[classSelect]-.2,speciesMeans[classSelect],weights[classSelect]+.2,speciesMeans[classSelect],col=speciesCols[classSelect],lwd=2)
  thisIntercepts<-sims[['intercepts']][,ii]
  thisBetas<-sims[['betas']][,ii]
  abline(mean(thisIntercepts),mean(thisBetas),lty=2)
}
dev.off()

dat2<-read.csv('tableScott.validation.csv',row.names=1,stringsAsFactors=FALSE)
#muegge<-dat2[dat2$study=='muegge',]
#pdf('test.pdf');
#plot(muegge$log.weight,muegge$log.otus)
#text(muegge$log.weight,muegge$log.otus,muegge$common)
#dev.off()
#head(muegge)
#dat2<-dat2[!is.na(dat2$log.otus),]
dat2$studySpecies<-sprintf('%s-%s',dat2$study,dat2$common)
speciesIds2<-1:length(unique(dat2$studySpecies))
names(speciesIds2)<-unique(dat2$studySpecies)
dat2$speciesId<-speciesIds2[dat2$studySpecies]
studyIds<-1:length(unique(dat2$study))
names(studyIds)<-unique(dat2$study)
dat2$classId<-studyIds[dat2$study]
dat2$group<-dat2$study
speciesIdClasses2<-sapply(names(speciesIds2),function(xx)dat2[dat2$studySpecies==xx,][1,'study'])
dat2$taxa<-tolower(apply(dat2[,c('phylum','class','order','family','genus','species')],1,paste,collapse=' '))

fit2<-stan(
  model_code=stanCode,
  data=list(
    otus=dat2$log.otus,
    n=nrow(dat2),
    ids=dat2$speciesId,
    nIds=max(dat2$speciesId),
    weights=tapply(dat2$log.weight,dat2$speciesId,mean),
    classes=tapply(dat2$classId,dat2$speciesId,'[[',1),
    nClasses=max(dat2$classId)
  ),
  control=list(adapt_delta=.99),
  chains=32,iter=100000,thin=100
)

plotFit<-function(fit,dat){
  sims<-extract(fit)
  speciesMeans<-apply(sims[['mus']],2,mean)
  weights<-tapply(dat$log.weight,dat$speciesId,mean)
  speciesCols<-rainbow(max(dat$speciesId))
  xlim<-range(dat$log.weight)
  ylim<-range(dat$log.otus)
  #plotting
  plot(dat$log.weight,dat$log.otus,pch=21,bg=speciesCols[dat$speciesId],col=rainbow(max(dat$classId))[dat$classId],lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',xlim=xlim,ylim=ylim)
  segments(weights-.2,speciesMeans,weights+.2,speciesMeans,col=speciesCols,lwd=2)
  for(ii in 1:max(dat$classId)){
    withAs(dat=dat[dat$classId==ii,],plot(dat$log.weight,dat$log.otus,pch=21,bg=speciesCols[dat$speciesId],col=NA,lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',main=dat[1,'group'],xlim=xlim,ylim=ylim,cex=2))
    classSelect<-sapply(1:max(dat$speciesId),function(x)dat[dat$speciesId==x,'classId'][1])==ii
    segments(weights[classSelect]-.2,speciesMeans[classSelect],weights[classSelect]+.2,speciesMeans[classSelect],col=speciesCols[classSelect],lwd=2)
    thisIntercepts<-sims[['intercepts']][,ii]
    thisBetas<-sims[['betas']][,ii]
    abline(mean(thisIntercepts),mean(thisBetas),lty=2)
  }
}

pdf('data2.pdf',width=12,height=12)
  plotFit(fit2,dat2)
dev.off()



stanCode<-' data{
  int<lower=0> n;
  real otus[n]; //log
  int<lower=1> ids[n];
  int<lower=0> nIds;
  vector[nIds] weights; //log
  int<lower=1> classes[nIds];
  int<lower=0> nClasses;
  int<lower=1> studies[n];
  int<lower=0> nStudies;
}
parameters {
  vector[nClasses] intercepts;
  //real<lower=0> metaSigma;
  vector<lower=0>[nClasses] classSigmas;
  vector<lower=0>[nStudies] sigmas;
  //real<lower=0> sigma;
  vector[nIds] rawMus;
 // vector[nClasses] rawBetas;
  real metaBeta;
  real<lower=0> metaBetaSd;
  vector[nStudies] rawStudyBetas;
  real<lower=0> metaStudyBetaSd;
  //real<lower=0> speciesVar;
}
transformed parameters {
  vector[n] mus;
  vector[nClasses] betas;
  vector[nStudies] studyBetas;
  vector[nIds] speciesMus;
  betas=metaBeta+metaBetaSd*rawBetas;
  studyBetas=metaStudyBetaSd*rawStudyBetas;
  //mus=intercepts[classes[ids]] + betas[classes[ids]] .* weights + rawMus[ids]*metaSigma;
  speciesMus=rawMus .* classSigmas[classes];
  mus = intercepts[classes[ids]] + studyBetas[studies]+betas[classes[ids]] .* weights[ids] + speciesMus[ids];
}
model {
  //speciesVar ~ gamma(1,3);
  //sigmas ~ gamma(1,speciesVar);
  classSigmas ~ gamma(1,.1);
  sigmas ~ gamma(1,.1);
  metaStudyBetaSd ~ gamma(1,.1);
  rawStudyBetas ~ normal(0,1);
  //sigma ~ gamma(1,.1);
  rawBetas ~ normal(0,1);
  metaBetaSd ~ gamma(1,.1);
  //metaSigma ~ gamma(1,3);
  rawMus ~ normal(0,1);
  //otus ~ normal(mus[ids],sigmas[classes[ids]]);
  otus ~ normal(mus,sigmas[studies]);
}
'
commonCols<-c('common','class','log.otus','log.weight','study','taxa')
dat3<-rbind(dat[,commonCols],dat2[,commonCols])
dat3$studySpecies<-sprintf('%s-%s',dat3$study,dat3$common)
speciesIds3<-1:length(unique(dat3$taxa))
names(speciesIds3)<-unique(dat3$taxa)
dat3$speciesId<-speciesIds3[dat3$taxa]
dat3$studyClass<-sprintf('%s-%s',dat3$study,dat3$class)
studyClassIds<-1:length(unique(dat3$studyClass))
names(studyClassIds)<-unique(dat3$studyClass)
#dat3$classId<-studyClassIds[dat3$studyClass]
dat3$group<-dat3$class
classIds<-1:length(unique(dat3$class))
names(classIds)<-unique(dat3$class)
dat3$classId<-classIds[dat3$class]
studyIds<-1:length(unique(dat3$study))
names(studyIds)<-unique(dat3$study)
dat3$studyId<-studyIds[dat3$study]
speciesIdClasses3<-sapply(names(speciesIds3),function(xx)dat3[dat3$studySpecies==xx,][1,'study'])

fit3<-stan(
  model_code=stanCode,
  data=list(
    otus=dat3$log.otus,
    n=nrow(dat3),
    ids=dat3$speciesId,
    nIds=max(dat3$speciesId),
    weights=tapply(dat3$log.weight,dat3$speciesId,mean),
    classes=tapply(dat3$classId,dat3$speciesId,'[[',1),
    nClasses=max(dat3$classId),
    studies=dat3$studyId,
    nStudies=max(dat3$studyId)
  ),
  control=list(adapt_delta=.99),
  chains=32,iter=5000,thin=5
)

pdf('trace3.pdf');print(traceplot(fit3,par=c('metaBeta','betas','metaSigma','sigma')));plot(fit3,par=c('metaBeta','betas'));dev.off()

plotFit2<-function(fit,dat){
  sims<-extract(fit)
  #speciesMeans<-apply(sims[['speciesMus']],2,mean)
  studyMeans<-apply(sims[['studyBetas']],2,mean)
  weights<-tapply(dat$log.weight,dat$speciesId,mean)
  speciesCols<-rainbow(max(dat$speciesId),alpha=.7)
  studyCols<-rainbow(max(dat$studyId),alpha=.7)
  xlim<-range(dat$log.weight)
  ylim<-range(dat$log.otus)
  dat$adjusted<-dat$log.otus-studyMeans[dat$studyId]
  #plotting
  plot(dat$log.weight,dat$log.otus,pch=21,bg=speciesCols[dat$speciesId],col=rainbow(max(dat$classId))[dat$classId],lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',xlim=xlim,ylim=ylim)
  plot(dat$log.weight,dat$adjusted,pch=21,bg=studyCols[dat$studyId],col=rainbow(max(dat$classId))[dat$classId],lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',xlim=xlim,ylim=ylim,main='Adjusted')
  #plot(weights,speciesMeans,col=speciesCols,ylab='Species mean OTUs',xlab='Weight')
  for(ii in 1:max(dat$classId)){
    withAs(dat=dat[dat$classId==ii,],plot(dat$log.weight,dat$log.otus,pch=21,bg=speciesCols[dat$speciesId],col=NA,lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',main=dat[1,'group'],xlim=xlim,ylim=ylim,cex=2))
    classSelect<-sapply(1:max(dat$speciesId),function(x)dat[dat$speciesId==x,'classId'][1])==ii
    #segments(weights[classSelect]-.2,speciesMeans[classSelect],weights[classSelect]+.2,speciesMeans[classSelect],col=speciesCols[classSelect],lwd=2)
    thisIntercepts<-sims[['intercepts']][,ii]
    thisBetas<-sims[['betas']][,ii]
    abline(mean(thisIntercepts),mean(thisBetas),lty=2)
    withAs(dat=dat[dat$classId==ii,],plot(dat$log.weight,dat$adjusted,pch=21,bg=studyCols[dat$studyId],col=NA,lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',main=sprintf('Adjusted %s',dat[1,'group']),xlim=xlim,ylim=ylim,cex=2))
  }
  betas<-sims[['studyBetas']]
  betas<-betas[sample(nrow(betas),1e3),]
  colnames(betas)<-sapply(1:max(dat3$studyId),function(xx)dat3[dat3$studyId==xx,'study'][1])
  vpPlot(rep(colnames(betas),each=nrow(betas)),unlist(betas),las=2,ylab='Study offset')
  sigmas<-sims[['sigmas']]
  sigmas<-sigmas[sample(nrow(betas),1e3),]
  colnames(sigmas)<-sapply(1:max(dat3$studyId),function(xx)dat3[dat3$studyId==xx,'study'][1])
  vpPlot(rep(colnames(sigmas),each=nrow(sigmas)),unlist(sigmas),las=2,ylab='Standard deviation')
  betas<-sims[['betas']]
  colnames(betas)<-sapply(1:max(dat3$classId),function(xx)dat3[dat3$classId==xx,'class'][1])
  betas<-cbind('meta'=sims[['metaBeta']],betas)
  betas<-betas[sample(nrow(betas),1e3),]
  vpPlot(factor(rep(colnames(betas),each=nrow(betas)),levels=colnames(betas)),unlist(betas),las=2,ylab='Slope')
  abline(h=0,lty=2,col='red')
}
pdf('data3.pdf',width=12,height=12)
  plotFit2(fit3,dat3)
dev.off()


if(FALSE){
  nSpecies<-20
  mus<-rnorm(nSpecies,3,1)
  weights<-rgamma(nSpecies,1,.1)
  beta<-.1
  ns<-replicate(nSpecies,sample(10,1))
  ids<-rep(1:length(ns),ns)
  idClasses<-c(1:5,sample(1:5,length(ns)-5,TRUE)) #make sure at least 1 each
  classes<-idClasses[ids]
  betas<-rnorm(length(unique(classes)),beta,.05)
  x<-mapply(function(mu,weight,beta,n)rnorm(n,mu+beta*weight,1),mus,weights,betas,ns,SIMPLIFY=FALSE)
  xWeights<-rep(weights,sapply(x,length))
  vpPlot(ids,unlist(x),offsetXArgs=list(width=.3),pch=21,bg=rainbow(max(ids))[ids])
  plot(xWeights,unlist(x),pch=21,bg=rainbow(max(ids))[ids])
  fit<-stan(
    model_code=stanCode,
    data=list(
      otus=unlist(x),
      n=length(unlist(x)),
      ids=ids,
      nIds=length(x),
      weights=weights,
      classes=idClasses,
      nClasses=length(unique(classes))
    ),
    control=list(adapt_delta=.95),
    chains=6,iter=3000 # get a few extra samples as an example
  )
  print(fit)

  plot(fit)
  print(traceplot(fit))
  sims<-extract(fit)
  hist(sims[['metaBeta']])
}
