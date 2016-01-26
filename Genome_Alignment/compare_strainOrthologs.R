#!/usr/bin/Rscript 

options(warn=1)

library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)
require(RColorBrewer)

original.parameters=par()
options(width=9999)

myDF <- read.csv('sequenceIdentity.txt', header=T, sep='\t') #print(names(myDF)) #[1] 'Ortholog_strain1'    'Ortholog_strain2'    'ortholog_similarity'
myDF <- myDF[1:3] #only interested in first 3 columns

means.df = ddply(myDF, .(Ortholog_strain1, Ortholog_strain2), summarize, averageIdentity=round(mean(ortholog_similarity, na.rm=TRUE), 2))
print(means.df)
minimum <- min(means.df$averageIdentity)

my_palette <- colorRampPalette(c('red', 'yellow', 'green'))(n = 299)
c_all <- melt(means.df)
print(c_all)


##PRINT HEATMAP OF STRAIN IDENTITY SIMILARITY:
#--------------------------------------------------------------------------
png(filename=paste('../images/seqIdentity.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(corpus <- qplot(x=Ortholog_strain1, y=Ortholog_strain2, data=c_all, fill=value, geom='tile') +
					geom_text(aes(label=value), color='black', size=4) + 
					scale_fill_gradient(limits=c(minimum-5, 100), low='gold', high='green4') +
					labs(title='Campylobacter Pair-wise Sequence Identity Comparison', x=NULL, y=NULL) +
					guides(fill = guide_legend(title = 'Sequence\nSimilarity %', title.theme = element_text(size=15, angle = 0))) + theme(legend.text=element_text(size=12))  +
					theme(axis.text.x=element_text(angle=45, size=14, hjust=1, colour='black'), axis.text.y=element_text(size=14, hjust=1, colour='black')) )
dev.off()
