#!/usr/bin/Rscript 

options(warn=1)

#library('data.table')
library(plyr)
library(ggplot2)
library(reshape2)
library(reshape)
library(scales)
require(utils)
library(dplyr)
require(RColorBrewer)

original.parameters=par()
options(width=9999)

myDF <- read.csv('sequenceIdentity.txt', header=T, sep='\t') 
#print(names(myDF)) #[1] "Ortholog_strain1"    "Ortholog_strain2"    "ortholog_similarity" "Region_strain1"      "Region_strain2"     
myDF <- myDF[1:3] #only interested in first 3 columns

means.df = ddply(myDF, .(Ortholog_strain1, Ortholog_strain2), summarize, averageIdentity=mean(ortholog_similarity, na.rm=TRUE) )
print(means.df)

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
c_all <- melt(means.df)
print(c_all)

#means.df = ddply(myDF, .(Ortholog_strain1, Ortholog_strain2), function(x) mean(x[,'ortholog_similarity']))
#print(means.df)
#c_df <- cast(melt(myDF), Ortholog_strain1 ~ Ortholog_strain2 | myDF$ortholog_similarity, c(mean))
#print(c_df)


#count(myDF, c('Ortholog_strain1', 'Ortholog_strain2'))

#cell_stats = ddply(
#    .data = myDF #use the ANT data
#    , .variables = .(Ortholog_strain1,Ortholog_strain2) #uses each combination of cue and flanker
#    , .fun = function(x){ #apply this function to each combin. of cue & flanker
#        to_return = data.frame(
#            , acc = mean(myDF$ortholog_similarity)
#        )
#        return(to_return)
#    }
#    , .progress = 'text'
#)


##PRINT GENE COUNT HISTOGRAM: Graph depicts contribution of strain's genes from plasmids (if applicable) and chromosomes separately. 
#--------------------------------------------------------------------------
#png(filename=paste('seqIdentity.png', sep=''), width=3750,height=2750,res=300)
#par(mar=c(9.5,4.3,4,2))
#print(corpus <- qplot(x=Phylum1, y=Phylum2, data=c_all, fill=value, geom='tile') +
#					geom_text(aes(label=Significance), color="black", size=6) + 
#					scale_fill_gradient2(limits=c(-1, 1), low='brown4', high='cadetblue4') +
#					labs(title='Pair-wise Spearman Coefficient for All Environments', x=NULL, y=NULL) +
#					guides(fill = guide_legend(title = 'Phylum', title.theme = element_text(size=15, angle = 0))) + theme(legend.text=element_text(size=12))  +
#					theme(axis.text.x=element_text(angle=45, size=14, hjust=1, colour='black'), axis.text.y=element_text(size=14, hjust=1, colour='black')) )
#dev.off()

#print(h <- ggplot(myDF, aes(x=Filename, stat='bin')) + geom_bar() +
#i		labs(title='Gene Count by Campy Strain DNA Source', x='Campy. Strain', y='Gene Count\n') +
#    	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=15), text = element_text(size=18)) +
#		theme(axis.text.x=element_text(angle=45, size=16, hjust=1), axis.text.y=element_text(size=16), legend.position='none', plot.title = element_text(size=22)) )
###########################################################################
#
###PRINT GENE COUNT HISTOGRAM: Graph depicts contribution of strain's genes from plasmids (if applicable) and chromosomes together 
##--------------------------------------------------------------------------
#png(filename=paste('images/Pangenome_Histogram.png', sep=''), width=3750,height=2750,res=300)
#par(mar=c(9.5,4.3,4,2))
#print(h <- ggplot(myDF, aes(x=Strain, stat='bin', fill=factor(Filename))) + geom_bar() +
#		labs(title='Gene Count by Strain Pangenome', x='Campylobacter Strains', y='Gene Count\n') +
#    	guides(title.theme = element_text(size=15, angle = 90)) + theme(legend.text=element_text(size=15), text = element_text(size=18)) +
#		theme(axis.text.x=element_text(angle=45, size=16, hjust=1), axis.text.y=element_text(size=16), legend.position='none', plot.title = element_text(size=22)) )
###########################################################################
#
#dev.off()
