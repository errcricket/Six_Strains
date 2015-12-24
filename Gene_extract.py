'''
Created:	E. Reichenberger
Date:	11.20.2015

Purpose: Extract a section of a sequence if given a start and stop position. Additionally, this script can also convert the extracted sequence an amino acid sequence and compare this to a referenced sequence if available. Lastly, if a match is not present, the script will look for a match nucleotide by nucleotide.

Input: Genome sequence (nucleotide), amino acid sequence from specific gene (along with start/stop position)
Output: uv_gene.fa

Working python version: Python 3.4.3 :: Anaconda 2.3.0 (64-bit) :: Biopython 1.6
'''

import Bio
import os
from Bio import GenBank
from Bio import SeqIO
from Bio import SeqFeature

#genome_file = 'Campy3194c/Campy3194c.gb'
genome_file = 'Campy1147c/Campy1147c.gb'

def match(sequence, index):
	print('\nmatch')
	print(index, sequence)

def get_aminoAcid(nucleotide_sequence, amino_sequence, start_position, stop_position, string=''):
	index = -1

	while index < 3:
		print(index)
		extracted_sequence = nucleotide_sequence[start_position+index:stop_position+index]
		stringF = str(extracted_sequence.translate(table=11))
		stringC = str(extracted_sequence.complement().translate(table=11))
		stringRC = str(extracted_sequence.reverse_complement().translate(table=11))

		if stringF == amino_sequence:
			match(stringF, index)
			break
		if stringC  == amino_sequence:
			match(stringC, index)
			break
		if stringRC  == amino_sequence:
			match(stringRC, index)
			break
			
		#uncomment to see actual translations
		print('\ntranslate:', stringF, '\n', 'complement:', stringC, '\n', 'reverse_complement:', stringRC, '\n')
		#print('Bacteria table translate:', stringBT, '\n')

		index+=1

record = next(SeqIO.parse(genome_file, 'genbank'))
sequence = record.seq

start  = 49163 
stop = 50149 
AA = 'MDNTQKYAAIDLKSFYASVECILRKLDPLNTNLVVADESRTEKTICLAVSPALRSYNISGRLRLFELIQKVKTINYERLKIAKYFSAKSYNHLELINNPNLELDYIVAKPRMSTYIDYSSKIYSIYLKYFDPKDIHIYSIDEVFIDLTPYIKHYKLSADKLIENILFEILKTTQITATAGIGTNLYLAKIAMDILAKKQNINKDGLCIGYLDEMLYRRKLWQHTPINDFWRIGKGYATKLKSIGINNMGDLARYSLNNEDKLYQIFGVNTELLIDHAWGFESCTMQAIKEYKSKHISKVMAKVLPKPYSFKKARNMLKEIVDHMVRAN'
get_aminoAcid(sequence, AA, start, stop)
