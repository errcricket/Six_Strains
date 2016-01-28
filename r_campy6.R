#!/usr/bin/Rscript 

library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)

original.parameters=par()
options(width=9999)
options(warn=1)

myDF <- read.csv('Output/campy6_corpus_cds.txt', header=T, sep='\t') 
print(names(myDF)) #[1] 'Filename'   'Strain'     'DNA_Source' 'Locus_Tag'  'Product'    'Transl_Tbl' 'Note'       'Seq_AA'     'Protein_ID'

myDF <- myDF[c(1, 3, 2, 5)]
myDF <- data.table(myDF) #had factors, not strings?

p <- strsplit(as.character(myDF$Product), split = ',') #split columns w/ multiple products into multiple rows
newDF <- data.frame(Filename=rep(myDF$Filename, sapply(p, length)), DNA_Source=rep(myDF$DNA_Source, sapply(p, length)), Strain = rep(myDF$Strain, sapply(p, length)), Product = unlist(p))
myDF <- data.table(newDF) #had factors, not strings?

#for some reason, I need to run this twice.
p <- strsplit(as.character(myDF$Product), split = ';') #split columns w/ multiple product_id into multiple rows
newDF <- data.frame(Filename=rep(myDF$Filename, sapply(p, length)), DNA_Source=rep(myDF$DNA_Source, sapply(p, length)), Strain = rep(myDF$Strain, sapply(p, length)), Product = unlist(p))
myDF <- data.table(newDF) #had factors, not strings?

##PRINT GENE COUNT HISTOGRAM: Graph depicts contribution of strain's genes from plasmids (if applicable) and chromosomes separately. 
png(filename=paste('images/Plasmid_Chromo_Histogram.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(h <- ggplot(myDF, aes(x=Filename, stat='bin')) + geom_bar() +
		labs(title='Gene Count by Campylobacteria Strain DNA Source', x='Campylobacter Strain', y='Gene Count\n') +
   	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=15), text = element_text(size=18)) +
		theme(axis.text.x=element_text(angle=45, size=16, hjust=1), axis.text.y=element_text(size=16), legend.position='none', plot.title = element_text(size=22)) )
##--------------------------------------------------------------------------

##PRINT GENE COUNT HISTOGRAM: Graph depicts contribution of strain's genes from plasmids (if applicable) and chromosomes together 
dt <- myDF[, .(countStrain = .N), by = c('Strain', 'Filename', 'DNA_Source')][order(Strain, Filename, DNA_Source)]
dt[, yval := cumsum(countStrain) - 0.5 * countStrain, by = Strain] # add the y-values for the plot

png(filename=paste('images/Pangenome_Histogram.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(z <- ggplot(dt, aes(x = Strain, y = countStrain, fill = Filename)) + 
		geom_bar(stat = 'identity') + 
		#geom_text(data = dt, angle=90, aes(x = Strain, y = yval, label = paste(DNA_Source, paste('\n[', countStrain, ']', sep=''),  sep=' '), size=2)) +
		geom_text(data = dt, angle=90, nudge_y=50, aes(x = Strain, y = yval, label = paste(DNA_Source, paste('\n[', countStrain, ']', sep=''),  sep=' '), label_size=2)) +
		labs(title='[Gene Count] by Campylobacter Strain (Pangenome)', x='Campylobacter Strains', y='Gene Count\n') +
    	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=15), text = element_text(size=18)) +
		theme(axis.text.x=element_text(angle=45, size=16, hjust=1), axis.text.y=element_text(size=16), legend.position='none', plot.title = element_text(size=22)) )
##--------------------------------------------------------------------------

##PRINT GENE COUNT HEATMAP: Graph depicts contribution of strain's genes from plasmids (if applicable) and chromosomes together 
gc <- myDF[, .(geneCount = .N), by = c('Strain', 'Product')][order(Strain, Product)]
c_all = melt(gc)

png(filename=paste('images/gene_heatMap.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(corpus <- qplot(x=Product, y=Strain, data=c_all, fill=value, geom='tile') +
		scale_fill_gradient(low='deepskyblue4', high='gold', limit=c(0,max(c_all$value)), name='Expression') +
		labs(title='Gene Comparison Across Campylobacter Strains(Pangenome)', x='Gene Products', y='Strains') +
		theme(axis.text.x=element_text(angle=90, size=5, hjust=1), axis.text.y=element_text(size=9),  plot.title = element_text(size=22)) )

dev.off()
##--------------------------------------------------------------------------

write.table(myDF, file='Output/productSplit_campy6_corpus_cds.txt', quote=F, sep='\t', row.names=F)
