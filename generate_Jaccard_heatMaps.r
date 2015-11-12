#!/usr/bin/Rscript # change the working directory in R to the place where you have the input file.  options(warn=1)

library('plyr')
library(ggplot2)
library(reshape2)

original.parameters=par()
options(width=9999)

myDF <- read.csv('Output/jaccard_r.txt', header=T, sep='\t')
myDF$Jaccard_C <- round(100*myDF$Jaccard_C, 2)

c_all = melt(myDF)
print(c_all)


png(filename=paste('images/Jaccard_Heatmap.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
#print(corpus <- qplot(x=Strain1, y=Strain2, data=c_all, fill=value, geom='tile', colour='black', name=NULL) +
print(corpus <- qplot(x=Strain1, y=Strain2, data=c_all, fill=value, geom='tile') +
               geom_text(aes(label=value), color='black', size=4) +
					scale_fill_gradient(low='cadetblue4', high='seagreen2', limit=c(80,100)) +
					theme(panel.grid.minor = element_line(colour = 'black')) + 
               labs(title='Pair-Wise Jaccard Similarity Coefficient\nAcross 7 Campylobacter Strains', x=NULL, y=NULL) +
               guides(fill = guide_legend(title = 'Jaccard\nSimilarity\nCoefficient (%)', title.theme = element_text(size=15, angle = 0))) + theme(legend.text=element_text(size=12))  +
               theme(axis.text.x=element_text(angle=45, size=14, hjust=1, colour='black'), axis.text.y=element_text(size=14, hjust=1, colour='black')) )
dev.off()
