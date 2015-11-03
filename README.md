# Six Campy Strains Project
---------------------------


###The following scripts are are used to extract CDS information from 7 Campylobacter strains.
---------------------------------------------------------------------------------------------


#####Step 1. python extractGB.py 
######Extract data from genebank files, output saved to campy6_corpus_cds.txt

#####Step 2. Rscript r_campy6.R  
######Script to create simple histogram of gene counts (by strain pangenome, by strain chromosome/plasmid)

#####Step 3. python gene_comparison.py 
######Script does cursory statistics, calculates Jaccard Similarity Coefficient (JSC) (uncommenting print statements & saving stdio to file called jaccard_r.txt, & find genes unique to a particular strain as well as comparing genes between two strains. Outputs statistics to campy_statitics.txt & unique genes to Unique_Genes.txt.

#####Step 4. Rscript generate_Jaccard_heatMaps.r 
######Creates a (very lackluster) heatmap of JSCs. Colors currently hard to decipher between values

#####Additional Notes: As it currently stands, the genbank files are downloaded to a directory (of the same name). For this specific project, the files were split up into multiple files if a plasmid was present, or if the genome has not been fully assembled.
