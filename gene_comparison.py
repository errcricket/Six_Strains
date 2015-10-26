'''
Author:	E. Reichenberger
Date:		10.26.2015

Purpose: Compare the campy strains for unique genes (Product), gene counts, ...(to be done w/ & w/out plasmids)
'''

import sys
import os
import Bio

def compute_jaccard_index(set_1, set_2):
    n = len(set_1.intersection(set_2))
    return n / float(len(set_1) + len(set_2) - n)
    #return round(n / float(len(set_1) + len(set_2) - n), 4)

campy_dic = {}
strain_dic = {}

#headers <- 'Filename"   "Strain"     "DNA_Source" "Locus_Tag"  "Product"    "Protein_ID" "Strand"     "Transl_Tbl" "Seq_AA'
with open('campy6_corpus_cds.txt', 'r') as inputFile:
	inputFile.readline()
	lines = inputFile.readlines()

	for line in lines:
		line = line.replace('\n', '')
		sLine = line.split('\t')

		if sLine[1] not in strain_dic:
			strain_dic[sLine[1]] = []
		strain_dic[sLine[1]].append(sLine[4]) #appends gene product

		if sLine[0] not in campy_dic:
			campy_dic[sLine[0]] = []
		campy_dic[sLine[0]].append(sLine[4]) #appends gene product
	
		
with open('campy_statistics.txt', 'w') as outputFile:
	outputFile.write('File\tTotal_Gene_Count\tTotal_Unique_Genes\n')
	for c in campy_dic:
		outputFile.write(c + '\t' + str(len(campy_dic[c])) + '\t' + str(len(set(campy_dic[c]))) + '\n')
		print(c, len(campy_dic[c]), len(set(campy_dic[c])))


#########Computing Jaccard Coefficient
with open('jaccard_r.txt', 'w') as outputFile:
	outputFile.write('Strain1\tStrain2\tJaccard_C\n')
	#Compute Jaccard Similarity Coefficient across the pangenomes for each strain
	print('set1\tset2\tJaccard Similarity Coefficient')
	with open('jaccard_campyPG.txt', 'w') as outputFile_1:
		for c in campy_dic:
			string = ''
			for r in campy_dic:
				jaccard = compute_jaccard_index(set(campy_dic[c]), set(campy_dic[r]))
			#	print(c, r, jaccard)
				string = string + '\t' + str(round(jaccard, 3))
				string = string.lstrip()
			outputFile_1.write(string + '\n')

	#Compute Jaccard Similarity Coefficient for the files (Plasmid/Chrom) are separated)
	with open('jaccard_campyStrain.txt', 'w') as outputFile_2:
		for c in strain_dic:
			string = ''
			for r in strain_dic:
				jaccard = compute_jaccard_index(set(strain_dic[c]), set(strain_dic[r]))
				print(c, r, jaccard)
				string = string.lstrip()
				string = string + '\t' + str(round(jaccard, 3))
				#string = string + '\t' + str(jaccard)
				outputFile.write(c + '\t' + r + '\t' + str(jaccard) + '\n')
			outputFile_2.write(string + '\n')



'''
Campy1147c/Campy1147c_Chrom.gbf

Campy1147q/Campy1147q_Chrom_1.gbf
Campy1147q/Campy1147q_Chrom_2.gbf
Campy1147q/Campy1147q_Chrom_3.gbf

Campy1188c/Campy1188c_Chrom.gbf
Campy1188c/Campy1188c_Plasmid.gbf

Campy1246c/Campy1246c_Chrom.gbf
Campy1246c/Campy1246c_Plasmid.gbf

Campy1285c/Campy1285c_Chrom.gbf

Campy14076c/Campy14076c_Chrom.gbf

Campy3194c/Campy3194c_Chrom.gbf
Campy3194c/Campy3194c_Plasmid.gbf
'''
