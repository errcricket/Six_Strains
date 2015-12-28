#!/usr/bin/Rscript 

options(warn=1)

library(plyr)
library(ggplot2)
library(reshape2)

original.parameters=par()
options(width=9999)

myDF <- read.csv('Campy3194c_analysisOutput.txt', header=T, sep='\t') 
print(names(myDF)) #[1] 'Accession' 'Locus_Tag' 'Gene'      'Identity.' 'RFP'       'index'    

means.df = ddply(myDF, .(Accession), summarize, averageIdentity=round(mean(Identity, na.rm=TRUE), 2)) #find average identity for each accession number
print(means.df)
low <- min(myDF$Identity)-1
low2 <- low - .15
width = length(myDF$Identity)/2

subA1 <- subset(myDF, Accession == 'RM3194_chrom')
min1 <- min(subA1$Identity)
subA2 <- subset(myDF, Accession == 'RM3194_plasmid')
min2 <- min(subA2$Identity)

text <- paste('Average identity: ', means.df$Accession[1], ': ', means.df$averageIdentity[1], '% ', means.df$Accession[2], ': ', means.df$averageIdentity[2], '%', sep='')
text2 <- paste('Lowest identity: ', means.df$Accession[1], ': ', min1, '% ', means.df$Accession[2], ': ', min2, '%', sep='')

png(filename=paste('images/CDSIdentity.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(corpus <- ggplot(myDF, aes(Locus_Tag, Identity)) +
					geom_point(aes(colour=Accession), size=3) + 
					labs(title='CDS Identity Comparison', x='Locus Tag', y='Identity Percentage', size=8) +
					guides(fill = guide_legend(title = 'Accession Number', title.theme = element_text(size=15, angle = 0))) + 
					theme(legend.text=element_text(size=12))  + 
					annotate('text', label=text, x=width, y=low, size=5, colour='lightsteelblue4') +
					annotate('text', label=text2, x=width, y=low2, size=5, colour='lightsteelblue4') +
					theme(axis.text.x=element_blank(), axis.text.y=element_text(size=16, hjust=1, colour='black')) )
dev.off()
