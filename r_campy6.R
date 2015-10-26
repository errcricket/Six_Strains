#!/usr/bin/Rscript 

options(warn=1)

#library('data.table')
library('plyr')
library(ggplot2)
library(reshape2)
library(scales)

original.parameters=par()
options(width=9999)


myDF <- read.csv('campy6_corpus_cds.txt', header=T, sep='\t') 
print(names(myDF))
print(length(names(myDF)))

#[1] "FileName"   "Locus_Tag"  "Product"    "Protein_ID" "Strand"     "Transl_Tbl" "Seq_AA"    

png(filename=paste('Histogram.png', sep=''), width=750,height=550,res=72)
par(mar=c(9.5,4.3,4,2))
print(h <- ggplot(myDF, aes(x=FileName, stat='bin')) + geom_bar() +
		labs(title='Gene Count by Campy Strain', x='Campy. Strain', y='Gene Count\n') +
    	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=13)) +
		theme(axis.text.x=element_text(angle=90, size=14, hjust=1)) )
    	#guides(fill = guide_legend(title = 'Amino Acid', title.theme = element_text(size=15, angle = 0))) + theme(legend.text=element_text(size=13))  +
#    theme_minimal() + theme(text=element_text(size=14)) + theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
#   theme(panel.grid.major = element_line(size = .3, color = 'grey'), axis.line = element_line(size=.5, color = 'black'), text = element_text(size=16)) +
#	scale_y_continuous(label=scientific_format()) )
#, stat='bin', fill=factor(Locus_Tag)) 
#+ geom_histogram() +# facet_wrap(~ variable) + #scale_fill_brewer(palette='Blues') +


uLocus = sort(unique(unlist(myDF$Locus_Tag, use.names = FALSE)))
png(filename=paste('XHistogram.png', sep=''), width=750,height=550,res=72)
par(mar=c(9.5,4.3,4,2))
print(x <- qplot(factor(FileName), data=myDF, geom='bar', fill=factor(length(uLocus))))
#print(x <- qplot(factor(FileName), data=myDF, geom='bar', fill=factor(Locus_Tag)))
dev.off()
