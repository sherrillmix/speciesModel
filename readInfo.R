info<-read.csv('data/islandGut discovery metadata - map.tsv.csv',row.names=1,stringsAsFactors=FALSE)
info$weight<-info$Weight_to_use
info2<-read.csv('data/islandGut validation metadata - map.tsv.csv',row.names=1,stringsAsFactors=FALSE)
