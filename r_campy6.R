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
#print(names(myDF)) #[1] "Filename"   "Strain"     "DNA_Source" "Locus_Tag"  "Product"    "Transl_Tbl" "Note"       "Seq_AA"     "Protein_ID"

##PRINT GENE COUNT HISTOGRAM: Graph depicts contribution of strain's genes from plasmids (if applicable) and chromosomes separately. 
#--------------------------------------------------------------------------
png(filename=paste('Plasmid_Chromo_Histogram.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(h <- ggplot(myDF, aes(x=Filename, stat='bin')) + geom_bar() +
		labs(title='Gene Count by Campy Strain DNA Source', x='Campy. Strain', y='Gene Count\n') +
    	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=13)) +
		theme(axis.text.x=element_text(angle=90, size=14, hjust=1)) )
##########################################################################

##PRINT GENE COUNT HISTOGRAM: Graph depicts contribution of strain's genes from plasmids (if applicable) and chromosomes together 
#--------------------------------------------------------------------------
png(filename=paste('Pangenome_Histogram.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(h <- ggplot(myDF, aes(x=Strain, stat='bin', fill=factor(Filename))) + geom_bar() +
		labs(title='Gene Count by Strain Pangenome', x='Campy. Strain', y='Gene Count\n') +
    	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=13)) +
		theme(axis.text.x=element_text(angle=90, size=14, hjust=1)) )
##########################################################################

png(filename=paste('X_Pangenome_Histogram.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(h <- ggplot(myDF, aes(x=Strain, fill=Filename), label=Filename) + geom_bar() +
		labs(title='Gene Count by Strain Pangenome', x='Campy. Strain', y='Gene Count\n') +
		geom_text(aes(y=Strain, label=Filename), size=5) + 
    	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=13)) +
		theme(axis.text.x=element_text(angle=90, size=14, hjust=1), legend.position='none') )

png(filename=paste('Z_Pangenome_Histogram.png', sep=''), width=3750,height=2750,res=300)
par(mar=c(9.5,4.3,4,2))
print(p <- ggplot(myDF, aes(x = Strain, stat='bin', fill=Filename)) ) #+ geom_text(label=myDF$Filename))
#print(p <- ggplot(myDF, aes(x = Strain, y = Frequency)) +
#     geom_bar(aes(fill = Filename), stat="identity") +
 #    geom_text(aes(label = Filename, y=Strain), size = 3))



#png(filename=paste('X_Pangenome_Histogram.png', sep=''), width=3750,height=2750,res=300)
#par(mar=c(9.5,4.3,4,2))
#print(h <- ggplot(myDF, aes(x=Strain, fill=Filename)) + geom_bar(stat="identity", ymin=0, aes(y=value, ymax=value), position="dodge") +
#geom_text(aes(x=Strain, y=value, ymax=value, label=value, hjust=ifelse(sign(value)>0, 1, 0)), position = position_dodge(width=1)) + scale_y_continuous(labels = percent_format()))
#print(h <- geom_point(myDF, aes(x=Strain, stat='bin', colour=factor(Filename)), size = 3) + 
#		theme(axis.text.x=element_text(angle=90, size=5, hjust=1), legend.position='none') + 
#		labs(title='Gene Count by Strain Pangenome', x='Campy. Strain', y='Gene Count\n') + 
#		guides(col = guide_legend(ncol = 5)) + theme(legend.text=element_text(size=13)) +
#		theme(axis.text.x=element_text(angle=90, size=14, hjust=1)) )
	

#	geom_text(aes(label=ENVIRONMENT), size=5, angle=90, hjust=1.3) +  theme(text = element_text(size=14)) + theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
#		theme(axis.text.y=element_text(size=13, colour='black'))




    	#guides(fill = guide_legend(title = 'Amino Acid', title.theme = element_text(size=15, angle = 0))) + theme(legend.text=element_text(size=13))  +
#    theme_minimal() + theme(text=element_text(size=14)) + theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
#   theme(panel.grid.major = element_line(size = .3, color = 'grey'), axis.line = element_line(size=.5, color = 'black'), text = element_text(size=16)) +
#	scale_y_continuous(label=scientific_format()) )
#, stat='bin', fill=factor(Locus_Tag)) 
#+ geom_histogram() +# facet_wrap(~ variable) + #scale_fill_brewer(palette='Blues') +


#uLocus = sort(unique(unlist(myDF$Locus_Tag, use.names = FALSE)))
#png(filename=paste('XHistogram.png', sep=''), width=3750,height=2750,res=300)
#par(mar=c(9.5,4.3,4,2))
#print(x <- qplot(factor(Strain), data=myDF, geom='bar', fill=factor(length(uLocus))))
#print(x <- qplot(factor(Filename), data=myDF, geom='bar', fill=factor(Locus_Tag)))
dev.off()
