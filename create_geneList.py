'''
Created:	E. Reichenberger
Date:	3.7.2016

Purpose: To extact all gene_product names from GenBank files and place the information to a new text file. Output file should have the following format: Gene_product\tFeature_type\tCategory\n. #NOTE: Category will be an empty column.
'''

import Bio
import os
import sys
from Bio import GenBank
from Bio import SeqIO

#######################################################DEFINITIONS################################################################
def strip_it(string_name):
	stripers = ['[', ']', '\'', '"']
	for s in stripers:
		string_name = string_name.replace(s, '')
		string_name = string_name.lstrip()
	return string_name
#---------------------------------------------------------------------------------------------------------------------------------


######################################################### File List ##############################################################
fileList = []

with open('fileList.txt', 'r') as inputFile:
	for i in inputFile.readlines():
		i = i.replace('\n', '')
		fileList.append(i)

geneCorpus = {} #will contain list of text files (same as fileList, but .txt extension)
#---------------------------------------------------------------------------------------------------------------------------------

############################Searchable Genbank Files & write qualifying features to output file####################################
for f in fileList:
	string = ''
	with open(f, 'r') as handle:
		for record in SeqIO.parse(handle, 'genbank'):
			for feature in record.features:
				feature_type = feature.type
				product_name = strip_it(str(feature.qualifiers.get('product')))  #.lower() change case (will be needed later to find unique gene products)
				if product_name not in geneCorpus:
					geneCorpus[product_name] = {}
				geneCorpus[product_name] = feature_type
#---------------------------------------------------------------------------------------------------------------------------------


############################Create Corpus File ####################################
with open('Output/gene_corpus.txt', 'w') as outputFile:
	outputFile.write('\t'.join(['Gene_product', 'Feature_type', 'Category']) + '\n')
	for g in geneCorpus:
		string = '\t'.join([g, geneCorpus[g]]) + '\n'
		outputFile.write(string)
#---------------------------------------------------------------------------------------------------------------------------------
