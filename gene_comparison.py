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

campyFile_dic = {}
strain_dic = {}

#Create dictionary for each strain (including plasmids), & for each file (plas-chrom are separate)
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

		if sLine[0] not in campyFile_dic:
			campyFile_dic[sLine[0]] = []
		campyFile_dic[sLine[0]].append(sLine[4]) #appends gene product
	
		
#Calculate total & unique number of genes (for both strains and files)
with open('campy_statistics.txt', 'w') as outputFile:
	outputFile.write('File\tTotal_Gene_Count\tTotal_Unique_Genes\n')
	for c in campyFile_dic: #Chrome & Plasmids separate
		outputFile.write(c + '\t' + str(len(campyFile_dic[c])) + '\t' + str(len(set(campyFile_dic[c]))) + '\n')
		print(c, len(campyFile_dic[c]), len(set(campyFile_dic[c])))

	print('')
	for s in strain_dic: #Strains separate
		outputFile.write(s + '\t' + str(len(strain_dic[s])) + '\t' + str(len(set(strain_dic[s]))) + '\n')
		print(s, len(strain_dic[s]), len(set(strain_dic[s])))
	print('')
	print('')

	#initialize all dictionaries to hold unique entries (not the most efficient method)
	dic_genes = {}
	for c in campyFile_dic:
		if c == 'Campy1147q':
			dic_genes['Campy1147q'] = {'Campy1285c':[], 'Campy1188c':[], 'Campy1246c':[],'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
		if c == 'Campy1285':
			dic_genes['Campy1285c'] = {'Campy1188c':[], 'Campy1246c':[],'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
		if c == 'Campy1188c':
			dic_genes['Campy1188c'] = {'Campy1246c':[],'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
		if c == 'Campy1246':
			dic_genes['Campy1246c'] = {'Campy3194c':[],'Campy1147c':[],'Campy14076c':[]}
		if c == 'Campy3194c':
			dic_genes['Campy3194c'] = {'Campy1147c':[],'Campy14076c':[]}
		if c == 'Campy1147c':
			dic_genes['Campy1147c'] = {'Campy14076c':[]}


	#Find genes unique to each strain/file.
	strain_order = ['Campy1147q', 'Campy1285c', 'Campy1188c', 'Campy1246c', 'Campy3194c', 'Campy1147c', 'Campy14076c']
	dic_genes = {}

	outputFile.write('\n\n----------------------Unique Genes----------------------\n')
#	for index, o in enumerate(strain_order):
#	for index, o in enumerate(campyFile_dic):
#		index = 0 #activating this line will make index start from 0 each time. Deactivating will access only remaining items.
#		if o not in dic_genes:
#			while index < len(strain_order):
#				dic_genes[o][strain_order[index]] = []

#			print(index, o)
#			index+=1

	for index, c in enumerate(strain_dic):
		index2 = index+1
		while index2 < len(strain_order):
			for i in set(strain_dic[c]):
				if i not in set(strain_dic[strain_order[index2]]):
					dic_genes[c][strain_order[index2]].append(i)
			index2+=1
			print(dic_genes)
#		for i in set(campyFile_dic[c]):
#			for a in set(campyFile_dic):
#				if i not in set(campyFile_dic[a]) and i not in 
#					print(i)
		

'''
#########Computing Jaccard Coefficient
with open('jaccard_r.txt', 'w') as outputFile:
	outputFile.write('Strain1\tStrain2\tJaccard_C\n')
	#Compute Jaccard Similarity Coefficient across the pangenomes for each strain
	print('set1\tset2\tJaccard Similarity Coefficient')
	with open('jaccard_campyPG.txt', 'w') as outputFile_1:
		for c in campyFile_dic:
			string = ''
			for r in campyFile_dic:
				jaccard = compute_jaccard_index(set(campyFile_dic[c]), set(campyFile_dic[r]))
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