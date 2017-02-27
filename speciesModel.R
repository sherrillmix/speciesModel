library(rstan)
library(vipor)
#add a regression parameter while retaining groups
stanCode<-'
data{
  int<lower=0> n;
  real x[n]; 
  int<lower=1> ids[n];
  int<lower=0> nIds;
  vector<lower=0>[nIds] weights;
  int<lower=1> classes[nIds];
  int<lower=0> nClasses;
}
parameters {
  real metaMu;
  real<lower=0> metaSigma;
  real<lower=0> sigma;
  vector[nIds] rawMus;
  vector[nClasses] rawBetas;
  real metaBeta;
  real<lower=0> metaBetaSd;
}
transformed parameters {
  vector[nIds] mus;
  vector[nIds] betas;
  betas=metaBeta+metaBetaSd*rawBetas[classes];
  mus = metaMu + betas .* weights + rawMus*metaSigma;
}
model {
  rawBetas ~ normal(0,1);
  metaBetaSd ~ gamma(1,1);
  metaSigma ~ gamma(1,1);
  rawMus ~ normal(0,1);
  x ~ normal(mus[ids],sigma);
}
'
mus<-rnorm(10,3,1)
weights<-rgamma(10,1,.1)
beta<-.1
x<-mapply(function(mu,weight)rnorm(sample(10,1),mu+beta*weight,1),mus,weights,SIMPLIFY=FALSE)
xWeights<-rep(weights,sapply(x,length))
ids<-rep(1:length(x),sapply(x,length))
idClasses<-c(1:5,sample(1:5,length(x)-5,TRUE)) #make sure at least 1 each
classes<-idClasses[ids]
vpPlot(ids,unlist(x))
plot(xWeights,unlist(x))
fit<-stan(
  model_code=stanCode,
  data=list(
    x=unlist(x),
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
traceplot(fit)
traceplot(fit)
sims<-extract(fit)
hist(sims[['mu']])


