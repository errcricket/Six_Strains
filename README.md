# Six Campy Strains Project
---------------------------


###The following scripts are are used to extract CDS information from 7 Campylobacter strains.
---------------------------------------------------------------------------------------------


#####Step 1. python extractGB.py 
######Extract data from genebank files, output saved to campy6_corpus_cds.txt (Filename Strain DNA_Source Locus_Tag Product Protein_ID Strand Transl_Tbl Seq_AA)

#####Step 2. Rscript r_campy6.R  
######Script to create simple histogram of gene counts (by strain pangenome, by strain chromosome/plasmid)

#####Step 3. python gene_comparison.py 
######Script does cursory statistics, calculates Jaccard Similarity Coefficient (JSC) (uncommenting print statements & saving stdio to file called jaccard_r.txt), & finds genes unique to a particular strain as well as comparing genes between two strains. Outputs statistics to campy_statitics.txt & unique genes to Unique_Genes.txt, also creates text file containing extracted CDS info for each strain genbank file.

#####Step 4. Rscript generate_Jaccard_heatMaps.r 
######Creates a (very lackluster) heatmap of JSCs. Colors currently hard to decipher between values

#####Step 5. python extract_gene.py 
######This script, using start/stop position of gene within genome will find extract the nucleotide sequence from the ORIGIN sequence, convert it to a protein sequence, & compare that to the printed protein sequence of the gene (in CDS section)

#####Additional Notes: Script only seems to work with .gbf extensions. As it currently stands, the genbank files are downloaded to a directory (of the same name). For this specific project, the files were split up into multiple files if a plasmid was present, or if the genome had not been fully assembled. All text files will be placed in the "Output" directory -- should add code to check for directory, and if not present, creates it. Running Python 2.7.10 :: Anaconda 2.3.0 (64-bit).


