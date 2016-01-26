''' Author:	E. Reichenberger Date:		11.6.2015

Purpose: Compare Mauve alignment (post formatting with fasta_formatter.py) of 7 campy strains. Compute pair-wise sequence identity for all strains. This script does not consider a separate scoring scheme for gaps.

file organized as below where each grouping represents the alignment to one ortholog:
>0:1-1323:Campy1147c_10 +
atgaatccaagccaaatact
...
>15:1-1323:Campy1285c_10 +
atgaatcc

>0:1-323:Campy1147c_10 +
---aatccaagccaaatact
...
>15:11140:Campy1147q_12350 -
atgaatcc
'''

import os

##FIND SEQUENCE NAME: Some strains are lacking ortholog. If lacking, header == >0
##########################################################################
def return_strain_name(header):
	strains = ['Campy1147c', 'Campy1147q', 'Campy1188c', 'Campy1246c', 'Campy1285c', 'Campy14076c',  'Campy3194c', 'CP006702', 'CP007179', 'CP007181', 'PSU1', 'PSU15', 'PSU29', 'PSU31', 'PSU32', 'RM1529']
	header_strains = ['>0', '>1', '>2', '>3', '>4', '>5', '>6', '>7', '>8', '>9', '>10', '>11', '>12', '>13', '>14', '>15']
	index = header_strains.index(header)
	
	return strains[index]
#--------------------------------------------------------------------------

##CALCULATE SEQUENCE SIMILARITY: Find % shared nucleotides
##########################################################################
def calculate_sequence_similarity(dictionary, f_name):
	identity = 0

	for i in dictionary:
		for j in dictionary:
			track = 0
		
			for index, k in enumerate(dictionary[j]):
				if k == dictionary[i][index]:
					track+=1

			identity = 100*round(track/len(dictionary[i]), 4)

			outputString = i + '\t' + j + '\t' + str(round(identity, 4)) + '\n' #+ region_strain1 + '\t' + region_strain2 + '\n'
			f_name.write(outputString)
#--------------------------------------------------------------------------

dic = {}

##Populate dictionary
##########################################################################
with open('/home/cricket/Projects/Campy_6Strains/Genome_Alignment/Campy16_formatted_alignments', 'r') as inputFile:
	o_name = '/home/cricket/Projects/Campy_6Strains/Genome_Alignment/sequenceIdentity.txt'

	try: #file in append mode, need to delete each time script is run
		 os.remove(o_name)
	except OSError:
		 pass

	with open(o_name, 'a') as outputFile:
		outputFile.write('Ortholog_strain1\tOrtholog_strain2\tortholog_similarity\n')

		lines = inputFile.readlines()
		for index, line in enumerate(lines):
			line = line.replace('\n', '')
			
			if line.startswith('>'): #Get sequences into dictionary w/ header as key
				line = line.rsplit(':')[0]
				strain_name = return_strain_name(line)
				#print(strain_name)
				
				if strain_name not in dic:
					dic[strain_name] = ''
				dic[strain_name] = lines[index+1].replace('\n', '').upper()
			
			if len(dic) == 16: #Have all orthologs from each strain, now compare them.
				calculate_sequence_similarity(dic, outputFile)
				dic = {} #recycle dictionary
#--------------------------------------------------------------------------
