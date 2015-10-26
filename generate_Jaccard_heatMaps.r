#!/usr/bin/Rscript # change the working directory in R to the place where you have the input file.  options(warn=1)

library('plyr')
library(ggplot2)
library(reshape2)

original.parameters=par()
options(width=9999)

myDF <- read.csv('jaccard_r.txt', header=T, sep='\t')
print(myDF)

c_all = melt(myDF)
#png(filename=paste('Jaccard_HeatMap.png', sep=''), width=3750,height=2750,res=300)
#par(mar=c(9.5,4.3,4,2))
#print(corpus <- qplot(x=Strain1, y=Strain2, data=c_all, fill=value, geom='tile') +
#		geom_text(aes(label=Strain1), color='black', size=6) + 
#		#scale_colour_gradient(limits=c(0.83, 1), low='brown4', high='cadetblue4') +
#		#scale_colour_gradient2(limits=c(0.83, 1)) + #, low='brown4', high='cadetblue4') +
#		scale_fill_gradient2(limits=c(0.83, 1), low='brown4', high='cadetblue4') +
#		labs(title='Jaccard Similarity Coefficient\nfor Campylobacter Pan-Genome', x=NULL, y=NULL) +
#		guides(fill = guide_legend(title = 'Jaccard Similarity\nCoefficient', title.theme = element_text(size=15, angle = 0))) + theme(legend.text=element_text(size=12))  +
#		theme(axis.text.x=element_text(angle=45, size=14, hjust=1, colour='black'), axis.text.y=element_text(size=14, hjust=1, colour='black')) )
#dev.off()


#png(filename=paste('Jaccard_HeatMap.png', sep=''), width=3750,height=2750,res=300)
#par(mar=c(9.5,4.3,4,2))
#print(corpus <- qplot(x=Strain1, y=Strain2, data=c_all, fill=value, geom='tile')  +
#			geom_text(aes(label=Strain1), color='black', size=3) + 
#			#scale_fill_gradient2(limits=c(0.83, 1), low='brown4', high='cadetblue4') +
#			scale_colour_gradient2(name=NULL, low = 'lightpink3', mid = 'cadetblue4', high = 'white', midpoint = 0.93, space='rgb') + 
#			#scale_fill_gradient2(limits=c(0.8, 1), high='lightpink3', low='cadetblue4') +
#			labs(title='Jaccard Similarity Coefficient\nfor Campylobacter Pan-Genome', x=NULL, y=NULL) +
#			guides(fill = guide_legend(title = 'Jaccard Similarity\nCoefficient', title.theme = element_text(size=15, angle = 0))) + theme(legend.text=element_text(size=12)) + 
#			theme(axis.text.x=element_text(angle=45, size=14, hjust=1, colour='black'), axis.text.y=element_text(size=14, hjust=1, colour='black')) )
#dev.off()


# Fill gradients work much the same way
png(filename=paste('puke.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(p <- qplot(letters[1:5], 1:5, fill= c(-3, 3, 5, 2, -2)) + geom_bar() + scale_fill_gradient2('fill'))
dev.off()
#
#print(p <- ggplot(myDF, aes(as.character(Filename), stat='bin', fill=factor(Product))) + geom_histogram() + #facet_wrap(~ variable) + #scale_fill_brewe
#      labs(title='Gene Count by Campy Strain', x='Campylobacter Strain', y='Gene Count\n') +
#      guides(fill = guide_legend(title = 'Amino Acid', title.theme = element_text(size=15, angle = 0))) + theme(legend.text=element_text(size=13))  +
#      theme_minimal() + theme(text=element_text(size=14)) + theme(axis.text.x=element_text(angle=90, size=14), axis.text.y=element_text(size=14)) +
#      theme(panel.grid.major = element_line(size = .3, color = 'grey'), axis.line = element_line(size=.5, color = 'black'), text = element_text(size=16
#      theme(legend.position='none') + scale_y_continuous(label=scientific_format()) )

