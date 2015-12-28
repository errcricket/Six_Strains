# Six Campy Strains Project
---------------------------


###The following scripts are are used to extract CDS information from 7 Campylobacter strains.
---------------------------------------------------------------------------------------------


#####Step 1. python extractGB.py 
######Extract data (Filename, Strain, DNA_Source, Locus_Tag, Gene Product, Protein_ID, Strand, Transl_Tbl, Seq_AA) from genbank files. NOTE: Multiple gene products (separated by comma) exist for some CDS. This script does not account for this and is dealt with in the r_campy6.R script.
######**Input**: Genbank Files 
######**Output**: creates text file containing extracted CDS info for each genbank file and campy6_corpus_cds.txt 

#####Step 2. Rscript r_campy6.R  
######Script to create simple histograms of gene counts (by strain pangenome, by strain chromosome/plasmid) and a gene expression heatmap. This script separates gene products (separated by comma). 
######**Input**: campy6_corpus_cds.txt
######**Output**: productSplit_campy6_corpus_cds.txt, Plasmid_Chromo_Histogram.png, Pangenome_Histogram.png, gene_heatMap.png 

#####Step 3. python gene_comparison.py 
######Script does cursory statistics, calculates Jaccard Similarity Coefficient (JSC) & finds genes unique to a particular strain as well as comparing genes between two strains. 
######**Input**: productSplit_campy6_corpus_cds.txt 
######**Output**: campy_statistics.txt, jaccard_statistics.txt, Unique_Genes.txt 

#####Step 4. Rscript generate_Jaccard_heatMaps.r 
######Creates a heatmap of JSCs. 
######**Input**: jaccard_statistics.txt 
######**Output**: Jaccard_Heatmap.png 

#####Step 5. python extract_gene.py 
######This script, using start/stop position of gene within genome will find extract the nucleotide sequence from the ORIGIN sequence, convert it to a protein sequence, & compare that to the printed protein sequence of the gene (in CDS section) 
######**Input**: Genome sequence (nucleotide), amino acid sequence from specific gene (along with start/stop position) 
######**Output**: uv_gene.fa

#####Step 6. Align genomes in Mauve using default parameters
######**Input**: Genebank files (*.gb extension (does not seem to work with *.gbf) 
######**Output**: Export alignment file 

#####Step 7. python format_fasta.py
######This script, places individual reads on single line. Example:
######>header_information
######*aagagagacgtatcgatcgatcgat*
######*gatcgtactactatctacgtgtgatcgt*

######becomes:

######>header_information 
######*aagagagacgtatcgatcgatcgatgatcgtactactatctacgtgtgatcgt*
######**Input**: Alignment file 
######**Output**: Fasta file 

#####Step 8. python calculate_orthologIdentity.py
######This script, uses the formatted (step 7) output from Mauve ortholog alignment (across all 7 strains) to calculate sequence similarity (based on positional nucleotides). 
######**Input**: Formatted fasta file (alignment) 
######**Output**: sequenceIdentity.txt 

#####Step 9. Rscript compare_strainOrthologs.R
######This script, uses output from step 8 to plot heatmap of strain similarity.
######**Input**: sequenceIdentity.txt 
######**Output**: seqIdentity.png  

#####Step 10. iTol
######Create Phylo tree using iTol (http://itol.embl.de/)
######**Input**: Alignment tree file (from Mauve) 
######**Output**: Campy_tree.png

#####Step 11. analyse_annotation.py
######This script will open/parse genbank files, extract CDS features, and compare (e.g., calculate sequence identity) the amino acid translation to the sliced & translated nucleotide sequence. If a match or close match, the calculated identity, locus_tag, and accession number are written to *_analysisOutput.txt. If a match is not close, the nucleotide slice is written to a fasta file (*_unpairedAnnotations.fasta) where the header is comprised of '> ' + locus_tag. Should add argument calls to script in future.
######**Input**: Genbank file(s)
######**Output**: *_analysisOutput.txt, *unpairedAnnotations.fasta 

#####Step 12. annotation_accuracy.R
######Creates graph of annotation identities to access accuracy of gb file. Should add automated method to label graph in future. The graph may work better as a boxplot.
######**Input**: *_analysisOutput.txt
######**Output**: CDSIdentity.png

#####Additional Notes: Python script only seems to work with .gbf extensions (but Mauve needs .gb extension). As it currently stands, the genbank files are downloaded to a directory (of the same name). For this specific project, the files were split up into multiple files if a plasmid was present, or if the genome had not been fully assembled. All text files will be placed in the "**Output**:" directory -- should add code to check for directory, and if not present, creates it. Running Python 2.7.10 : Anaconda 2.3.0 (64-bit) : Biopython 1.6. 
#####


