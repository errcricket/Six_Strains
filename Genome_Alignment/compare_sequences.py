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

##CALCULATE SEQUENCE SIMILARITY: Find % shared nucleotides
##########################################################################
def calculate_sequence_similarity(dictionary):
	with open('/home/cricket/Projects/Campy_6Strains/Genome_Alignment/sequenceIdentity.txt', 'w') as outputFile:
		outputFile.write('Ortholog_strain1\tOrtholog_strain2\tortholog_similarity\tRegion_strain1\tRegion_strain2\n')

		for i in dictionary:
			for j in dictionary:
				track = 0
				for index, k in enumerate(dictionary[j]):
					if k == dictionary[i][index]:
						track+=1
				outputString = i + '\t' + j + '\t' + str(100*round(track/len(dictionary[i]), 4)) + '\n'
				outputFile.write(outputString)
	#			print(i, j, str(track/len(dictionary[i])))

dic = {}

with open('/home/cricket/Projects/Campy_6Strains/Genome_Alignment/small_Formatted_Ortholog_Campy7.alignments', 'r') as inputFile:

	lines = inputFile.readlines()
	for index, line in enumerate(lines):
		line = line.replace('\n', '')
		
		if line.startswith('>'): #Get sequences into dictionary w/ header as key
			strain = line.rsplit('_')[0].rsplit(':')[2]
			if strain not in dic:
				dic[strain] = ''
			dic[strain] = lines[index+1].replace('\n', '')
		
		if len(dic) == 7: #Have all 7 orthologs from each strain, now compare them.
			oString = ''
			calculate_sequence_similarity(dic)
			dic = {} #recycle dictionary

