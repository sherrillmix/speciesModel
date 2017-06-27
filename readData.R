library(dnar)
nCut<-1000
nReads<-1:10


otu<-read.table('data/otu_table.txt.gz',sep='\t',comment='',skip=1,header=TRUE,stringsAsFactors=FALSE,row.names=1)
info<-read.csv('data/islandGut discovery metadata - map.tsv.csv',row.names=1,stringsAsFactors=FALSE)
taxa<-otu[,ncol(otu)]
otu<-otu[,-ncol(otu)]
n<-apply(otu,2,sum)
#prop<-apply(otu,2,function(xx)xx/sum(xx))
isEnough<-n>nCut
rareN<-do.call(cbind,lapply(nReads,function(xx){apply(otu[,isEnough],2,rareEquation,nCut/2,xx)}))
colnames(rareN)<-nReads
write.csv(rareN,'work/rareN.csv')

speciesAbund<-apply(otu,2,function(xx)sort(xx[xx>0]))
save(speciesAbund,file='work/speciesAbund.Rdat')


otu2<-read.table('data/otu_table_validation.txt.gz',sep='\t',comment='',skip=1,header=TRUE,stringsAsFactors=FALSE,row.names=1)
info2<-read.csv('data/islandGut validation metadata - map.tsv.csv',row.names=1,stringsAsFactors=FALSE)
taxa2<-otu2[,ncol(otu2)]
otu2<-otu2[,-ncol(otu2)]
#prop2<-apply(otu2,2,function(xx)xx/sum(xx))
n2<-apply(otu2,2,sum)
isEnough2<-n2>nCut
rareN2<-do.call(cbind,lapply(1:10,function(xx){apply(otu2[,isEnough2],2,rareEquation,nCut/2,xx)}))
colnames(rareN2)<-1:10
write.csv(rareN2,'work/rareN2.csv')

speciesAbund2<-apply(otu2,2,function(xx)sort(xx[xx>0]))
save(speciesAbund2,file='work/speciesAbund2.Rdat')

pdf('out/rare.pdf')
plot(1,1,type='n',xlim=c(1,ncol(rareN)),ylim=range(rareN),log='y',xlab='Counts required',ylab='Number of species (2500 reads)',las=1)
apply(rareN,1,function(xx)lines(1:ncol(rareN),xx,col='#00000044'))
dev.off()
