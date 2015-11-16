#!/usr/bin/Rscript 

options(warn=1)

#library('data.table')
library('plyr')
library(ggplot2)
library(reshape2)
library(scales)
library(data.table)

original.parameters=par()
options(width=9999)

myDF <- read.csv('Output/campy6_corpus_cds.txt', header=T, sep='\t') 
#print(names(myDF)) #[1] 'Filename'   'Strain'     'DNA_Source' 'Locus_Tag'  'Product'    'Transl_Tbl' 'Note'       'Seq_AA'     'Protein_ID'

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
myDF <- myDF[c(1, 2, 3)]
myDF <- data.table(myDF) #had factors, not strings?

dt <- myDF[, .(countStrain = .N), by = c('Strain', 'Filename', 'DNA_Source')][order(Strain, Filename, DNA_Source)]
#dt <- myDF[, .(countStrain = .N), by = c('Strain', 'DNA_Source')][order(Strain, DNA_Source)]

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
