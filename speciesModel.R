library(rstan)
library(vipor)

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
  real metaMu;
  real<lower=0> metaSigma;
  vector<lower=0>[nIds] sigmas;
  vector[nIds] rawMus;
  vector[nClasses] rawBetas;
  real metaBeta;
  real<lower=0> metaBetaSd;
  real <lower=0> speciesVar;
}
transformed parameters {
  vector[nIds] mus;
  vector[nClasses] betas;
  betas=metaBeta+metaBetaSd*rawBetas;
  mus = metaMu + betas[classes] .* weights + rawMus*metaSigma;
}
model {
  sigmas ~ gamma(1,speciesVar);
  rawBetas ~ normal(0,1);
  metaBetaSd ~ gamma(1,2);
  metaSigma ~ gamma(1,2);
  rawMus ~ normal(0,1);
  otus ~ normal(mus[ids],sigmas[ids]);
}
'



dat<-read.csv('data2.csv',row.names=1)
dat$speciesId<-as.numeric(as.factor(dat$host_species_common_name))
dat$classId<-as.numeric(as.factor(dat$class))

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
  chains=8,iter=3000,thin=10
)
print(fit)

plot(fit)
print(traceplot(fit))
sims<-extract(fit)
hist(sims[['metaBeta']])

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
