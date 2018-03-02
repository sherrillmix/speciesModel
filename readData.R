library(dnar)
nCut<-200
nReads<-1:10


otu<-read.table('data/otu_table.txt.gz',sep='\t',comment='',skip=1,header=TRUE,stringsAsFactors=FALSE,row.names=1)
info<-read.csv('data/islandGut discovery metadata - map.tsv.csv',row.names=1,stringsAsFactors=FALSE)
info$weight<-info$Weight_to_use
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

deblur<-read.table('data/deblur.txt.gz',sep='\t',comment='',skip=1,header=TRUE,stringsAsFactors=FALSE,row.names=1)
deblurAbund<-apply(deblur[,colnames(deblur)!='Consensus.Lineage'],2,function(xx)sort(xx[xx>0]))
deblur2<-read.table('data/deblur2.txt.gz',sep='\t',comment='',skip=1,header=TRUE,stringsAsFactors=FALSE,row.names=1)
deblurAbund2<-apply(deblur2[,colnames(deblur2)!='Consensus.Lineage'],2,function(xx)sort(xx[xx>0]))
deblur4<-read.table('data/deblur4.txt.gz',sep='\t',comment='',skip=1,header=TRUE,stringsAsFactors=FALSE,row.names=1)
deblurAbund4<-apply(deblur4[,colnames(deblur4)!='Consensus.Lineage'],2,function(xx)sort(xx[xx>0]))
save(deblurAbund,deblurAbund2,deblurAbund4,file='work/deblurAbund.Rdat')

dada<-read.table('data/dada2.otu.tsv.gz',sep='\t',comment='',skip=0,header=TRUE,stringsAsFactors=FALSE)
dadaAbund<-apply(dada,2,function(xx)sort(xx[xx>0]))
save(dadaAbund,file='work/dadaAbund.Rdat')
meta<-read.csv('islandGut discovery metadata - map.tsv.csv',stringsAsFactors=FALSE)
rownames(meta)<-meta$X.SampleID
goodIds<-meta$X.SampleID[meta$toKeep=='keep']
nDada<-apply(dada,2,sum)
#dadaIsEnough<-nDada>1000
dadaIsEnough<-names(nDada) %in% goodIds
dadaRareN<-apply(dada[,dadaIsEnough],2,rareEquation,500,1)
dadaWeights<-rbind(info[,'weight',drop=FALSE],info2[,'weight',drop=FALSE])[colnames(dada)[dadaIsEnough],'weight']
dadaSpecies<-sub('unknown_or_none $',' sp.',sub('acheta domestica','acheta domesticus',sub('sus scrofa domesticus','sus scrofa ',apply(rbind(info[,c('genus','species','subspecies')],info2[,c('genus','species','subspecies')])[colnames(dada)[dadaIsEnough],],1,paste,collapse=' '))))
dadaClass<-rbind(info[,c('class'),drop=FALSE],info2[,c('class'),drop=FALSE])[colnames(dada)[dadaIsEnough],]
dadaSpecies[grepl('^unknown_or_none',dadaSpecies)]<-NA
dadaCommon<-meta[colnames(dada)[dadaIsEnough],'common']
dadaDat<-data.frame(
  'weight'=log10(tapply(dadaWeights,dadaSpecies,mean)),
  'otu'=log10(tapply(dadaRareN,dadaSpecies,mean)),
  'class'=tapply(dadaClass,dadaSpecies,unique),
  'common'=tapply(dadaCommon,dadaSpecies,unique),
  stringsAsFactors=FALSE
)
summary(lm(otu~weight,dadaDat))

pdf('out/dadaSpeciesArea.pdf')
  plot(10^dadaDat$weight,10^dadaDat$otu,log='xy',xaxt='n',yaxt='n',xlab='Animal weight (g)',ylab='Number of OTUs (rarefied to 500 reads)',cex=.5)
  logAxis(1)
  logAxis(2,las=1)
  text(10^dadaDat$weight,10^dadaDat$otu,rownames(dadaDat),cex=.5)
  thisLm<-lm(otu~weight,dadaDat)
  fakeWeights<-seq(min(dadaDat$weight,na.rm=TRUE)-1,max(dadaDat$weight)+1,length.out=1000)
  preds<-predict(thisLm,data.frame('weight'=fakeWeights),interval='conf')
  lines(10^fakeWeights,10^preds[,'fit'])
  polygon(10^c(fakeWeights,rev(fakeWeights)),10^c(preds[,'lwr'],rev(preds[,'upr'])),border=NA,col='#00000022')
dev.off()


library(ggrepel)
library(scales)
scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
#classColors<-structure(
  #c("#FF0000FF", "#FF9900FF", "#CCFF00FF", "#33FF00FF", "#00FF66FF", "#00FFFFFF", "#0066FFFF", "#3300FFFF", "#CC00FFFF", "#FF0099FF"),
  #.Names = c("Actinopteri", "Aves", "Chondrichthyes", "Diplopoda", "Holothuroidea", "Insecta", "Malacostraca", "Mammalia", "Polychaeta", "Reptilia")
#)
classColors<-structure(
  c('#8dd3c7','#ffed6f','#bc80bd','#fccde5','#80b1d3','#fdb462','#b3de69','#fb8072','#bebada','#ccebc5'),
  .Names = c("Actinopteri", "Aves", "Chondrichthyes", "Diplopoda", "Holothuroidea", "Insecta", "Malacostraca", "Mammalia", "Polychaeta", "Reptilia")
)
pdf('out/dadaSpeciesArea.pdf',width=7,height=7,useDingbats=FALSE)
print(
  ggplot(dadaDat, aes(10^weight, 10^otu, label = sub(' ','\n',common))) +
    geom_smooth(method=lm, color='#00000033')+
    geom_text_repel(size=2.5,box.padding=.12,point.padding=.3,lineheight=.7,min.segment.length=.1,max.iter=3e4,nudge_y=.02,nudge_x=.1,color='#00000099') +
    #scale_x_continuous(trans='log10',label=scientific_10)+
    ylab('Number of sequence variants (rarefied to 500 reads)')+
    xlab('Animal weight (g)')+
    theme_classic(base_size = 16)+
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x,n=3),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(sides = "lb",short = unit(1,"mm"), mid = unit(1,"mm"), long = unit(2,"mm"))+
    geom_point(fill = classColors[dadaDat$class],pch=21,size=4)
)
dev.off()

