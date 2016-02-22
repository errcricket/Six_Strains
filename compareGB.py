'''
Created:	E. Reichenberger
Date:	2.22.2015

Purpose: To extact information from GenBank file and place the information to a new file. Output file should have the following format:
'Feature', 'Start', 'Stop', 'Product_name' (sep by tabs)

Output text files can be imported into Excel for side by side comparison of the annotations 
'''

import Bio
import os
import sys
from Bio import GenBank
from Bio import SeqIO

#######################################################DEFINITIONS################################################################
def strip_it(string_name):
	stripers = ['[', ']', '\'']
	for s in stripers:
		string_name = string_name.replace(s, '')
		string_name = string_name.lstrip()
	return string_name
#---------------------------------------------------------------------------------------------------------------------------------

###############################################File List & Searchable Qualifiers###################################################
fileList = []

with open('fileList.txt', 'r') as inputFile:
	for i in inputFile.readlines():
		i = i.replace('\n', '')
		fileList.append(i)

textList = [] #will contain list of text files (same as fileList, but .txt extension)
#---------------------------------------------------------------------------------------------------------------------------------

############################Searchable Genbank Files & write qualifying features to output file####################################
for f in fileList:
	extension = os.path.splitext(f)[1]
	nameOut = ''
	file_name = ''

	nameOut = f.replace(extension, '.txt')
	textList.append(nameOut) #add to list for 
	print(f, nameOut)

	string = ''
	with open(f, 'r') as handle:
		with open(nameOut, 'w') as outputFile:
			outputFile.write('\t'.join(['Feature', 'Start', 'Stop', 'Product']) + '\n')
			for record in SeqIO.parse(handle, 'genbank'):
				#print(len(record.features))

				for feature in record.features:
					if feature.type != 'gene':
						feature_type = feature.type

						feature_position = feature.location
						feature_position = strip_it(str(feature_position))
						position = feature_position.split('(')[0]

						start = position.split(':')[0]	
						stop = position.split(':')[1]	

						product_name = feature.qualifiers.get('product')
						product_name = strip_it(str(product_name))
						outputFile.write('\t'.join([feature_type, start, stop, product_name]) + '\n')
						print(feature_type, start, stop, product_name)
						print('')
#---------------------------------------------------------------------------------------------------------------------------------
