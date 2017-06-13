otu<-read.table('otu_table.txt.gz',sep='\t',comment='',skip=1)
n<-apply(otu,2,sum)
prop<-apply(otu,2,function(xx)xx/sum(xx))
