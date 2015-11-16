#!/usr/bin/Rscript 

options(warn=1)

library(ggplot2)
library(plyr)
library(data.table)

original.parameters=par()
options(width=9999)

myDF <- read.csv('Output/campy6_corpus_cds.txt', header=T, sep='\t') 
#print(names(myDF)) #[1] 'Filename'   'Strain'     'DNA_Source' 'Locus_Tag'  'Product'    'Transl_Tbl' 'Note'       'Seq_AA'     'Protein_ID'

myDF <- myDF[c(1, 3, 2, 5)]
print(names(myDF))
myDF <- data.table(myDF) #had factors, not strings?

p <- strsplit(as.character(myDF$Product), split = ',') #split columns w/ multiple products into multiple rows
newDF <- data.frame(Filename=rep(myDF$Filename, sapply(p, length)), DNA_Source=rep(myDF$DNA_Source, sapply(p, length)), Strain = rep(myDF$Strain, sapply(p, length)), Product = unlist(p))
myDF <- data.table(newDF) #had factors, not strings?

trim.leading <- function (x)  sub('^\\s+', '', x) #remove leading spaces from cell
myDF$Product <- trim.leading(myDF$Product)

##PRINT GENE COUNT HISTOGRAM: Graph depicts contribution of strain's genes from plasmids (if applicable) and chromosomes separately. 
#--------------------------------------------------------------------------
png(filename=paste('images/Plasmid_Chromo_Histogram.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(h <- ggplot(myDF, aes(x=Filename, stat='bin')) + geom_bar() +
		labs(title='Gene Count by Campy Strain DNA Source', x='Campy. Strain', y='Gene Count\n') +
    	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=15), text = element_text(size=18)) +
		theme(axis.text.x=element_text(angle=45, size=16, hjust=1), axis.text.y=element_text(size=16), legend.position='none', plot.title = element_text(size=22)) )
##########################################################################

##PRINT GENE COUNT HISTOGRAM: Graph depicts contribution of strain's genes from plasmids (if applicable) and chromosomes together 
#--------------------------------------------------------------------------

dt <- myDF[, .(countStrain = .N), by = c('Strain', 'Filename', 'DNA_Source')][order(Strain, Filename, DNA_Source)]
dt[, yval := cumsum(countStrain) - 0.5 * countStrain, by = Strain] # add the y-values for the plot

png(filename=paste('images/Pangenome_Histogram.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(z <- ggplot(dt, aes(x = Strain, y = countStrain, fill = Filename)) + 
		geom_bar(stat = 'identity') + 
		geom_text(data = dt, aes(x = Strain, y = yval, label = paste(DNA_Source, paste('[', countStrain, ']', sep=''),  sep=' '), size=5)) +
		labs(title='[Gene Count] by Campylobacter Strain (Pangenome)', x='Campylobacter Strains', y='Gene Count\n') +
    	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=15), text = element_text(size=18)) +
		theme(axis.text.x=element_text(angle=45, size=16, hjust=1), axis.text.y=element_text(size=16), legend.position='none', plot.title = element_text(size=22)) )

dev.off()
##########################################################################


write.table(myDF, file='Output/productSplit_campy6_corpus_cds.txt', quote=F, sep='\t', row.names=F)
