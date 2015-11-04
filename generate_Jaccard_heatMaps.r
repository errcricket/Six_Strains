#!/usr/bin/Rscript # change the working directory in R to the place where you have the input file.  options(warn=1)

library('plyr')
library(ggplot2)
library(reshape2)

original.parameters=par()
options(width=9999)

myDF <- read.csv('Output/jaccard_r.txt', header=T, sep='\t')
print(myDF)

breaks <- cut(min(myDF$Jaccard_C):max(myDF$Jaccard_C), 3)
c_all = melt(myDF)
nh_human = melt(myDF)


png(filename=paste('images/Puke.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(c <- ggplot(c_all, aes(Strain1, Strain2))+geom_tile(data=c_all, aes(fill=value), color='white')+
	scale_fill_gradient2(low='blue', high='red', mid='white', midpoint=0, limit=c(0.8,1),name='Jaccard\nSimilarity\nCoefficiet')+
	theme(axis.text.x = element_text(angle=45, vjust=1, size=11, hjust=1))+coord_equal())

png(filename='images/pukeII.png', width=700,height=550,res=72)
print(ip <- qplot(x=Strain1, y=Strain2, data=nh_human, fill=value, geom='tile') +
					#geom_text(aes(label=Significance), color='black', size=6) + 
					scale_fill_gradient2(limits=c(.8, 1), low='brown4', high='cadetblue4') +
					labs(title='Pair-wise Spearman Coefficient for Non-Human Gut Environments', x=NULL, y=NULL) +
					guides(fill = guide_legend(title = 'Phylum', title.theme = element_text(size=15, angle = 0))) + theme(legend.text=element_text(size=12))  +
					theme(axis.text.x=element_text(angle=45, size=14, hjust=1, colour='black'), axis.text.y=element_text(size=14, hjust=1, colour='black')) )
dev.off()
