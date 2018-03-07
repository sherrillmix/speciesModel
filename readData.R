library(dnar)
nCut<-200
if(!dir.exists('work'))dir.create('work')

source('readInfo.R')

otu<-read.table('data/otu_table.txt.gz',sep='\t',comment='',skip=1,header=TRUE,stringsAsFactors=FALSE,row.names=1)
taxa<-otu[,ncol(otu)]
otu<-otu[,-ncol(otu)]
n<-apply(otu,2,sum)
isEnough<-n>nCut
rareN<-apply(otu[,isEnough],2,rareEquation,nCut/2)
write.csv(cbind(names(rareN),'rare'=rareN),'work/rareN.csv',row.names=FALSE)

speciesAbund<-apply(otu,2,function(xx)sort(xx[xx>0]))
save(speciesAbund,file='work/speciesAbund.Rdat')


otu2<-read.table('data/otu_table_validation.txt.gz',sep='\t',comment='',skip=1,header=TRUE,stringsAsFactors=FALSE,row.names=1)
taxa2<-otu2[,ncol(otu2)]
otu2<-otu2[,-ncol(otu2)]
n2<-apply(otu2,2,sum)
isEnough2<-n2>nCut
rareN2<-apply(otu2[,isEnough2],2,rareEquation,nCut/2)
write.csv(cbind(names(rareN2),'rare'=rareN2),'work/rareN.csv',row.names=FALSE)

speciesAbund2<-apply(otu2,2,function(xx)sort(xx[xx>0]))
save(speciesAbund2,file='work/speciesAbund2.Rdat')


dada<-read.table('data/dada2.otu.tsv.gz',sep='\t',comment='',skip=0,header=TRUE,stringsAsFactors=FALSE)
dadaAbund<-apply(dada,2,function(xx)sort(xx[xx>0]))
save(dadaAbund,file='work/dadaAbund.Rdat')

goodIds<-rownames(info)[info$toKeep=='keep']
dadaRareN<-apply(dada[,colnames(dada) %in% goodIds],2,rareEquation,500,1)
write.csv(cbind(names(dadaRareN),'rare'=dadaRareN),'work/dadaRareN.csv',row.names=FALSE)

