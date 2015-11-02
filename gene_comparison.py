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
with open('campy6_corpus_cds.txt', 'r') as inputFile: inputFile.readline()
	lines = inputFile.readlines()

	for line in lines:
		line = line.replace('\n', '')
		sLine = line.split('\t')

		if sLine[1] not in strain_dic:
			strain_dic[sLine[1]] = []
		strain_dic[sLine[1]].append(sLine[4]) #appends gene product

		if sLine[0] not in campyFile_dic:
			campyFile_dic[sLine[0]] = []
		campyFile_dic[sLine[0]].append(sLine[4]) #appends gene product
#--------------------------------------------------------------------------
	
		
##CUSORY STATISTICS: Calculate total & unique number of genes (e.g., sets for both strains and files)
##########################################################################
#NOTE: uncomment '#print' for testing purposes
with open('campy_statistics.txt', 'w') as outputFile:
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
with open('campy_statistics.txt', 'a') as outputFile:
	outputFile.write('\n\nJACCARD SIMILARITY COEFFICIENT (by chromosome & plasmid)------------------\n')
#	outputFile.write('Strain1\tStrain2\tJaccard_C\n') #this is for r & will need to be in a separate file
	#print('set1\tset2\tJaccard Similarity Coefficient')

	orderString = ''
	for c in campyFile_dic: #Compute Jaccard Similarity Coefficient across the pangenomes for each strain
		string = c
		orderString = orderString + '\t' + c
		for r in campyFile_dic:
			jaccard = compute_jaccard_index(set(campyFile_dic[c]), set(campyFile_dic[r]))
			string = string + '\t' + str(round(jaccard, 3))
			#outputFile.write(c + '\t' + r + '\t' + str(jaccard) + '\n') #this is for r
		string = string.lstrip()
		outputFile.write(string + '\n')
	outputFile.write(orderString.lstrip() + '\n')

	outputFile.write('\n\nJACCARD SIMILARITY COEFFICIENT (by pangenome)------------------\n')
#	outputFile.write('Strain1\tStrain2\tJaccard_C\n') #this is for r & will need to be in a separate file

	orderString = ''
	for c in strain_dic: #Compute Jaccard Similarity Coefficient for the files (Plasmid/Chrom) are separated)
		string = c
		orderString = orderString + '\t' + c
		for r in strain_dic:
			jaccard = compute_jaccard_index(set(strain_dic[c]), set(strain_dic[r]))
			string = string + '\t' + str(round(jaccard, 3))
			#outputFile.write(c + '\t' + r + '\t' + str(jaccard) + '\n') #this is for r
		string = string.lstrip()
		outputFile.write(string + '\n')
	outputFile.write(orderString.lstrip() + '\n')
#--------------------------------------------------------------------------

##This function fills arrays (dictionary[primaryKey][s]) w/ genes that are found in strain_dic[primaryKey] but not in strain_dic[subkey]
def find_uniqueGenes(dictionary, primaryKey, strainList):
	dictionary[primaryKey] = {} #initialize dic w/ primary key
	for s in strainList:
		if s not in dictionary[primaryKey]:
			dictionary[primaryKey][s] = [] #initialize dic w/ primary key & subkey

		for i in set(strain_dic[primaryKey]):
			if i not in set(strain_dic[s]):
				dictionary[primaryKey][s].append(i)

	return dictionary
		
###################################################################################################
##FIND UNIQUE GENES: #Find genes unique to each strain/file.
strain_order = ['Campy1147q', 'Campy1285c', 'Campy1188c', 'Campy1246c', 'Campy3194c', 'Campy1147c', 'Campy14076c']

for index, s in enumerate(strain_order):
	unique_genes = {}
	unique_genes = find_uniqueGenes(unique_genes, s, strain_order[index+1:]) #blank after column indicates last entry
	print(unique_genes)
	print('')
	
#unique_genes = find_uniqueGenes(unique_genes, 'Campy1147q', strain_order[1:]) #blank after column indicates last entry
#['Campy1147q']
#['Campy1285c']
#['Campy1188c']
#['Campy1246c']
#['Campy3194c']
#['Campy1147c']


##########################################################################
#--------------------------------------------------------------------------
	#initialize all dictionaries to hold unique entries (not the most efficient method)
#	dic_genes = {}
#	for c in campyFile_dic:
#		if c == 'Campy1147q':
#			dic_1147q['Campy1147q'] = {'Campy1285c':[], 'Campy1188c':[], 'Campy1246c':[],'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
#		if c == 'Campy1285':
#			dic_1285c['Campy1285c'] = {'Campy1188c':[], 'Campy1246c':[],'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
#		if c == 'Campy1188c':
#			dic_1188c['Campy1188c'] = {'Campy1246c':[],'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
#		if c == 'Campy1246':
#			dic_1246c['Campy1246c'] = {'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
#		if c == 'Campy3194c':
#			dic_3194c['Campy3194c'] = {'Campy1147c':[],'Campy14076c':[]}
#		if c == 'Campy1147c':
#			dic_1147c['Campy1147c'] = {'Campy14076c':[]}

#dic_1147q['Campy1147q'] = {'Campy1285c':[], 'Campy1188c':[], 'Campy1246c':[],'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
#dic_1285c['Campy1285c'] = {'Campy1188c':[], 'Campy1246c':[],'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
#dic_1188c['Campy1188c'] = {'Campy1246c':[],'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
#dic_1246c['Campy1246c'] = {'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
#dic_3194c['Campy3194c'] = {'Campy1147c':[],'Campy14076c':[]}
#dic_1147c['Campy1147c'] = {'Campy14076c':[]}



	#Find genes unique to each strain/file.
#	strain_order = ['Campy1147q', 'Campy1285c', 'Campy1188c', 'Campy1246c', 'Campy3194c', 'Campy1147c', 'Campy14076c']
	#dic_genes = {}

#	outputFile.write('\n\n----------------------Unique Genes----------------------\n')

#dic_1147q = NestedDict()
#dic_1285c = NestedDict()
#dic_1188c = NestedDict()
#dic_1246c = NestedDict()
#dic_3194c = NestedDict()
#dic_1147c = NestedDict()
#
#reusable_dic = NestedDict()
#
##def find_uniqueGenes(dic, master_strain, slave_strain):
##	for s in strain_dic:
##		for t in strain_dic:
##			if s != t:
##				dic[s][t]
#
##for s in strain_dic:
#dic_1147q['Campy1285c']['Campy1188c']['Campy1246c']['Campy3194c']['Campy1147c']['Campy14076c']
#dic_1285c['Campy1188c']['Campy1246c']['Campy3194c']['Campy1147c']['Campy14076c']
#dic_1188c['Campy1246c']['Campy3194c']['Campy1147c']['Campy14076c']
#dic_1246c['Campy3194c']['Campy1147c']['Campy14076c']
#dic_3194c['Campy1147c']['Campy14076c']
#dic_1147c['Campy14076c']
#print(dic_1147q)
#print(dic_1147q[1])
#
##strain_order = ['Campy1147q', 'Campy1285c', 'Campy1188c', 'Campy1246c', 'Campy3194c', 'Campy1147c', 'Campy14076c']
##for s in strain_dic:
##	find_uniqueGenes(reusable_dic, 
##	resable_dic = NestedDict()
#
#	
#
#
#
##	for index, o in enumerate(strain_order):
##	for index, o in enumerate(campyFile_dic):
##		index = 0 #activating this line will make index start from 0 each time. Deactivating will access only remaining items.
##		if o not in dic_genes:
##			while index < len(strain_order):
##				dic_genes[o][strain_order[index]] = []
#
##			print(index, o)
##			index+=1
#
##	for index, c in enumerate(strain_dic):
##		index2 = index+1
##		while index2 < len(strain_order):
##			for i in set(strain_dic[c]):
##				if i not in set(strain_dic[strain_order[index2]]):
##					dic_genes[c][strain_order[index2]].append(i)
##			index2+=1
##			print(dic_genes)
##		for i in set(campyFile_dic[c]):
##			for a in set(campyFile_dic):
##				if i not in set(campyFile_dic[a]) and i not in 
##					print(i)
#		
#
##	gene_dic = defaultdict(lambda: defaultdict(dict))
##	for c in strain_dic:
##		if c not in gene_dic:
##			#gene_dic[c] = {}
##			for g in gene_dic:
##				if g != c:
##					for i in set(gene_dic[c]):	
##						if i not in gene_dic[g]:
##							gene_dic[c][g].append(i)
#
##print(gene_dic)
#
#
##missing = {}
##def find_uniqueGenes(strain1, strain2):
##	for strain1, strain2 in permutations(genes, 2): 
##		 missing.setdefault(strain1, {})[strain2] = genes[strain1] - genes[strain2]
##
##	  
###	strain_order = ['Campy1147q', 'Campy1285c', 'Campy1188c', 'Campy1246c', 'Campy3194c', 'Campy1147c', 'Campy14076c']
##genes = { 
##    'Campy1147q':  set(strain_dic['Campy1147q']),
##    'Campy1285c':  set(strain_dic['Campy1285c']),
##    'Campy1188c':  set(strain_dic['Campy1188c']),
##    'Campy1246c':  set(strain_dic['Campy1246c']),
##    'Campy3194c':  set(strain_dic['Campy3194c']),
##    'Campy1147c':  set(strain_dic['Campy1147c']),
##    'Campy14076c': set(strain_dic['Campy14076c']),
##}
##
##for c in genes:
##	for g in genes:
##		find_uniqueGenes(genes[c], genes[g])
##
##for m in missing:
##	print(m, len(missing[m]), len(m))
###print(missing)
###v
###for strain1, strain2 in permutations(genes, 2): 
###    missing.setdefault(strain1, {})[strain2] = genes[strain1] - genes[strain2]
###
###assert len(missing) == 3
###assert missing['strain1']['strain2'] == set([1, 2]) 
###assert missing['strain1']['strain3'] == set([1, 2, 3]) 
###assert missing['strain2']['strain1'] == set([4, 5]) 
###assert missing['strain2']['strain3'] == set([3, 4]) 
###assert missing['strain3']['strain2'] == set([6, 7]) 
###eggs = NestedDict()
###eggs[1][2][3][4][5]
###print(eggs)
##
##
##
##
