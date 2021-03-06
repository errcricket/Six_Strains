'''
Author:	E. Reichenberger
Date:		10.26.2015

Purpose: Compare the campy strains (file campy6_corpus_cds.txt) for unique genes (Product), gene counts, ...(to be done w/ & w/out plasmids)
Three files will be created. 
1. campy_statistics.txt contains very cursory information such as the name of the file, the total number of genes, & the total number of unique genes.
2. 
3. This file contains the Jaccard Similarity Coefficients for the strains and for the files. 

headers for campy6_corpus_cds.txt: 'Filename"   "Strain"     "DNA_Source" "Locus_Tag"  "Product"    "Protein_ID" "Strand"     "Transl_Tbl" "Seq_AA'
'''

import sys
import os
#import Bio
from collections import defaultdict
from itertools import permutations
import itertools

#######################################################DEFINITIONS################################################################
class NestedDict(dict):
	def __getitem__(self, key):
		if key in self: return self.get(key)
		return self.setdefault(key, NestedDict())

def compute_jaccard_index(set_1, set_2):
    n = len(set_1.intersection(set_2))
    return n / float(len(set_1) + len(set_2) - n)
    #return round(n / float(len(set_1) + len(set_2) - n), 4)
#--------------------------------------------------------------------------


##CREATE DICTIONARIES: Will use for all code, includes dictionary for each strain (including plasmids), & for each file (plas-chrom are separate)
##########################################################################
campyFile_dic = {} #holds genes for each chromosome & plasmids separately (even if from the same strain)
strain_dic = {} #holds genes for each strain (plasmid & chromosome collectively)
strains = [] #will keep a record of all strains (pangenome)
strain_files = [] 

with open('Output/productSplit_campy6_corpus_cds.txt', 'r') as inputFile: 
	inputFile.readline()
	lines = inputFile.readlines()

	for line in lines:
		line = line.replace('\n', '')
		sLine = line.split('\t')

		if sLine[2] not in strains:
			strains.append(sLine[2]) #PANGENOME array: corresponds to Strain column in productSplit_campy6_corpus_cds.txt

		if sLine[2] not in strain_dic: #PANGENOME: corresponds to Strain column in productSplit_campy6_corpus_cds.txt
			strain_dic[sLine[2]] = []
		strain_dic[sLine[2]].append(sLine[3]) #appends gene product

		if sLine[0] not in strain_files: 
			strain_files.append(sLine[0]) #CHROM/PLASMID array: corresponds to Filename column in productSplit_campy6_corpus_cds.txt

		if sLine[0] not in campyFile_dic: #CHROM/PLASMID: corresponds to Filename column in productSplit_campy6_corpus_cds.txt
			campyFile_dic[sLine[0]] = []
		campyFile_dic[sLine[0]].append(sLine[3]) #appends gene product
#--------------------------------------------------------------------------


##CUSORY STATISTICS: Calculate total & unique number of genes (e.g., sets for both strains and files) #NOTE: uncomment '#print' for testing purposes
##########################################################################
with open('Output/campy_statistics.txt', 'w') as outputFile:
	outputFile.write('CURSORY STATISTICS------------------\n')
	outputFile.write('File\tTotal_Gene_Count\tTotal_Unique_Genes\n')
	for c in campyFile_dic: #Chrome & Plasmids separate
		outputFile.write(c + '\t' + str(len(campyFile_dic[c])) + '\t' + str(len(set(campyFile_dic[c]))) + '\n')
		#print(c, len(campyFile_dic[c]), len(set(campyFile_dic[c])))

	#print('')
	for s in strain_dic: #Strains separate
		outputFile.write(s + '\t' + str(len(strain_dic[s])) + '\t' + str(len(set(strain_dic[s]))) + '\n')
		#print(s, len(strain_dic[s]), len(set(strain_dic[s])))
	#print('\n\n')
#--------------------------------------------------------------------------

##COMPUTING JACCARD COEFFICIENT: Compute Jaccard Similarity Coefficient across each strain pangenome & each entity (plasmid, chromosome)
###################################################################################################
with open('Output/campy_statistics.txt', 'a') as outputFile:
	outputFile.write('\n\nJACCARD SIMILARITY COEFFICIENT (by chromosome & plasmid)------------------\n')
	with open('Output/jaccard_statistics.txt', 'w') as J_outputFile:

		orderString = ''

		for c in strain_files: #want output in this order #Compute Jaccard Similarity Coefficient for the files (Plasmid/Chrom) are separated)
			string = c
			orderString = orderString + '\t' + c
			#print(orderString)
			for c2 in strain_files: #want output in this order #Compute Jaccard Similarity Coefficient for the files (Plasmid/Chrom) are separated)
				jaccard = compute_jaccard_index(set(campyFile_dic[c]), set(campyFile_dic[c2]))
				string = string + '\t' + str(round(jaccard, 3))

				#J_outputFile.write(c + '\t' + c2 + '\t' + str(jaccard) + '\n') #this is for r, but only want pangenome
			string = string.lstrip()
			outputFile.write(string + '\n')
		outputFile.write(orderString.lstrip() + '\n')

		outputFile.write('\n\nJACCARD SIMILARITY COEFFICIENT (by pangenome)------------------\n')
		J_outputFile.write('Strain1\tStrain2\tJaccard_C\n') #this is for r & will need to be in a separate file

		orderString = ''

		for c in strains: #want output in this order #Compute Jaccard Similarity Coefficient across the pangenomes for each strain
			string = c
			orderString = orderString + '\t' + c
			for c2 in strains: #want output in this order #Compute Jaccard Similarity Coefficient across the pangenomes for each strain
				print(c, c2)
				jaccard = compute_jaccard_index(set(strain_dic[c]), set(strain_dic[c2]))
				string = string + '\t' + str(round(jaccard, 3))
				J_outputFile.write(c + '\t' + c2 + '\t' + str(jaccard) + '\n') #this is for r
			string = string.lstrip()
			outputFile.write(string + '\n')
		outputFile.write(orderString.lstrip() + '\n')
#--------------------------------------------------------------------------


###################################################################################################
##FIND UNIQUE GENES: #Find genes unique to each strain/file.
##This function fills arrays (dictionary[primaryKey][s]) w/ genes that are found in strain_dic[primaryKey] but not in strain_dic[subkey]
def find_uniqueGenes(dictionary, primaryKey, strainList, completeStrainList):
	dictionary[primaryKey] = {} #initialize dic w/ primary key
	total_genes = []

	for s in completeStrainList: #get (almost all) gene products into single list
		if s != primaryKey:
			total_genes.append(strain_dic[s])

	for s in strainList: #get gene products found in strain_dic[primaryKey] that is not found in remaining strains
		if s not in dictionary[primaryKey]:
			dictionary[primaryKey][s] = [] #initialize dic w/ primary key & subkey

		for i in set(strain_dic[primaryKey]):
			if i not in set(strain_dic[s]):
				dictionary[primaryKey][s].append(i)

	# Look for unique genes across all strains save the strain in question
	dictionary[primaryKey]['Unique_to_all_strains'] = [] #initialize dic w/ primary key & subkey

	temp = list(itertools.chain(*total_genes)) #Cat gene entries from 6 (out of 7) strains, include unique entries only
	total_genes = list(set(temp))

	for i in set(strain_dic[primaryKey]):
		if i not in total_genes:
			dictionary[primaryKey]['Unique_to_all_strains'].append(i)
		
	uniqueList = dictionary[primaryKey]['Unique_to_all_strains']
	return dictionary, uniqueList
#--------------------------------------------------------------------------
		

###################################################################################################
##PRINT UNIQUE GENES: #Find genes unique to each strain/file.
with open('Output/Unique_Genes.txt', 'w') as outputFile:
	for index, s in enumerate(strains):
		outputFile.write('\n\n\subsubsection*{The following genes are unique to strain ' + s + '}\n')
		unique_genes = {}
		unique_genes, uList = find_uniqueGenes(unique_genes, s, strains[index+1:], strains) #blank after column indicates last entry
		outputFile.write(str(len(uList)) + ' '  + str(uList))
#--------------------------------------------------------------------------
