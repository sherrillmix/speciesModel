if(!file.exists('work/rareN.csv'))source('readData.R')
rare<-read.csv('work/rareN.csv',row.names=1)
rare2<-read.csv('work/rareN2.csv',row.names=1)
info<-read.csv('data/islandGut discovery metadata - map.tsv.csv',row.names=1,stringsAsFactors=FALSE)
info2<-read.csv('data/islandGut validation metadata - map.tsv.csv',row.names=1,stringsAsFactors=FALSE)
info$otu<-rare[rownames(info),'X1']
info2$otu<-rare2[rownames(info2),'X1']
info$study<-'bushman'
info$weight<-info$Weight_to_use
info$nRead<-info$filteredReadCount
info2$nRead<-info2$readCount
sharedCols<-intersect(colnames(info),colnames(info2))
combo<-rbind(info[,sharedCols],info2[info2$study!='bushman',sharedCols])
combo<-combo[combo$nRead>200&!is.na(combo$weight),]
combo$speciesName<-tolower(paste(combo$genus,combo$species))
combo$log.weight<-log10(combo$weight)
combo$log.otus<-log10(combo$otu)
weightSpeciesOrder<-tapply(combo$log.weight,combo$speciesName,mean)
combo$speciesId<-as.numeric(factor(combo$speciesName,levels=names(sort(weightSpeciesOrder))))
combo$classId<-as.numeric(as.factor(combo$class))
combo$studyId<-as.numeric(factor(combo$study,levels=c('bushman',unique(combo$study[combo$study!='bushman']))))


