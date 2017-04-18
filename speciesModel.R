library(rstan)
library(vipor)
library(dnar)

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
  //vector<lower=0>[nIds] sigmas;
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
  sigma ~ gamma(1,3);
  rawBetas ~ normal(0,1);
  metaBetaSd ~ gamma(1,3);
  metaSigma ~ gamma(1,3);
  rawMus ~ normal(0,1);
  otus ~ normal(mus[ids],sigma);
}
'



dat<-read.csv('tableScott.discovery.csv',row.names=1,stringsAsFactors=FALSE)
#CAREFUL excluding single sample classes for now
#dat<-dat[ave(dat$class,dat$class,FUN=length)>1,]
speciesIds<-1:length(unique(dat$host_species_common_name))
names(speciesIds)<-unique(dat$host_species_common_name)
dat$speciesId<-speciesIds[dat$host_species_common_name]
classIds<-1:length(unique(dat$class))
names(classIds)<-unique(dat$class)
dat$classId<-classIds[dat$class]
speciesIdClasses<-sapply(names(speciesIds),function(xx)dat[dat$host_species_common_name==xx,][1,'class'])


plot(dat$log.weight,dat$log.otu,pch=21,bg=rainbow(max(dat$speciesId))[dat$speciesId],col=rainbow(max(dat$classId))[dat$classId],lwd=2)
meanWeight<-tapply(dat$log.weight,dat$host_species_common_name,mean)
meanOtu<-tapply(dat$log.otu,dat$host_species_common_name,mean)
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
dat2<-dat2[!is.na(dat2$log.otus),]
speciesIds2<-1:length(unique(dat2$common))
names(speciesIds2)<-unique(dat2$common)
dat2$speciesId<-speciesIds2[dat2$common]
classIds2<-1:length(unique(dat2$study))
names(classIds2)<-unique(dat2$study)
dat2$classId<-classIds[dat2$study]
speciesIdClasses2<-sapply(names(speciesIds2),function(xx)dat2[dat2$common==xx,][1,'class'])

fit2<-stan(
  model_code=stanCode,
  data=list(
    otus=dat$log.otus,
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

sims2<-extract(fit2)
speciesMeans<-apply(sims2[['mus']],2,mean)
weights<-tapply(dat2$log.weight,dat2$speciesId,mean)
speciesCols<-rainbow(max(dat2$speciesId))
xlim<-range(dat2$log.weight)
ylim<-range(dat2$log.otus)
pdf('data2.pdf',width=12,height=12)
plot(dat2$log.weight,dat2$log.otus,pch=21,bg=speciesCols[dat2$speciesId],col=rainbow(max(dat2$classId))[dat2$classId],lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',xlim=xlim,ylim=ylim)
segments(weights-.2,speciesMeans,weights+.2,speciesMeans,col=speciesCols,lwd=2)
for(ii in 1:max(dat2$classId)){
  withAs(dat=dat2[dat2$classId==ii,],plot(dat$log.weight,dat$log.otus,pch=21,bg=speciesCols[dat$speciesId],col=NA,lwd=2,las=1,xlab='Log weight',ylab='Log rarefied species',main=names(classIds)[ii],xlim=xlim,ylim=ylim,cex=2))
  classSelect<-speciesIdClasses==names(classIds)[ii]
  segments(weights[classSelect]-.2,speciesMeans[classSelect],weights[classSelect]+.2,speciesMeans[classSelect],col=speciesCols[classSelect],lwd=2)
  thisIntercepts<-sims2[['intercepts']][,ii]
  thisBetas<-sims2[['betas']][,ii]
  abline(mean(thisIntercepts),mean(thisBetas),lty=2)
}
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
