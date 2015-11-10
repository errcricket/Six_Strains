''' Author:	E. Reichenberger Date:		11.6.2015

Purpose: Compare Mauve alignment (post formatting with fasta_formatter.py) of 7 campy strains. Compute pair-wise sequence identity for all strains. This script does not consider a separate scoring scheme for gaps.

file organized as below where each grouping represents the alignment to one ortholog:
>0:1-1323:Campy1147c_10 +
atgaatccaagccaaatact
...
>6:1-1323:Campy1285c_10 +
atgaatcc

>0:1-323:Campy1147c_10 +
---aatccaagccaaatact
...
>6:11140:Campy1147q_12350 -
atgaatcc
'''

import os

##FIND SEQUENCE NAME: Some strains are lacking ortholog. If lacking, header == >0
##########################################################################
def return_strain_name(header):
	strains = ['Campy1147c', 'Campy1147q', 'Campy1188c', 'Campy1246c', 'Campy1285c', 'Campy14076c',  'Campy3194c']
	header_strains = ['>0', '>1', '>2', '>3', '>4', '>5', '>6']
	index = header_strains.index(header)
	
	return strains[index]
#--------------------------------------------------------------------------

##CALCULATE SEQUENCE SIMILARITY: Find % shared nucleotides
##########################################################################
def calculate_sequence_similarity(dictionary, f_name):
		strain1 = ''
		strain2 = ''
		region_strain1 = ''
		region_strain2 = ''
		identity = 0

		for i in dictionary:
			if len(i) == 2: #if strain lacks ortholog, header is only '>X' where X = number from 0-6
				strain1 = return_strain_name(i)
				region_strain1 = 'NA'

			else:
				strain1= i.rsplit(':')[2]
				region_strain1 = i.rsplit(':')[1]

			for j in dictionary:
				if len(j) == 2:
					strain2 = return_strain_name(j)
					region_strain2 = 'NA'

				else:
					strain2= j.rsplit(':')[2]
					region_strain2 = j.rsplit(':')[1]

				track = 0

				if len(i) == 2 or len(j) == 2: #if header lacking, strain is lacking ortholog, identity = 0
					identity = 0
			
				else:
					while len(dictionary[j]) < len(dictionary[i]):
						dictionary[j] = dictionary[j] + '-'
					while len(dictionary[i]) < len(dictionary[j]):
						dictionary[i] = dictionary[i] + '-'

					for index, k in enumerate(dictionary[j]):
						if k == dictionary[i][index]:
							track+=1

					identity = 100*round(track/len(dictionary[i]), 4)

				if i == j:
					identity = 100.0
#				if identity != 0: #
				outputString = strain1 + '\t' + strain2 + '\t' + str(identity) + '\t' + region_strain1 + '\t' + region_strain2 + '\n'
				f_name.write(outputString)
#--------------------------------------------------------------------------

dic = {}

##CALCULATE SEQUENCE SIMILARITY: Find % shared nucleotides
##########################################################################
with open('/home/cricket/Projects/Campy_6Strains/Genome_Alignment/Formatted_Ortholog_Campy7.alignments', 'r') as inputFile:
	o_name = '/home/cricket/Projects/Campy_6Strains/Genome_Alignment/sequenceIdentity.txt'

	try: #file in append mode, need to delete each time script is run
		 os.remove(o_name)
	except OSError:
		 pass

	with open(o_name, 'a') as outputFile:
		outputFile.write('Ortholog_strain1\tOrtholog_strain2\tortholog_similarity\tRegion_strain1\tRegion_strain2\n')

		lines = inputFile.readlines()
		for index, line in enumerate(lines):
			line = line.replace('\n', '')
			
			if line.startswith('>'): #Get sequences into dictionary w/ header as key
				line = line.rsplit('_')[0]
				if line not in dic:
					dic[line] = ''
				dic[line] = lines[index+1].replace('\n', '')
			
			if len(dic) == 7: #Have all 7 orthologs from each strain, now compare them.
				calculate_sequence_similarity(dic, outputFile)
				dic = {} #recycle dictionary
#--------------------------------------------------------------------------
