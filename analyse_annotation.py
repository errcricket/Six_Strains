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
from Bio.Seq import Seq

#genome_file = 'sequence.gb'
genome_file = 'Campy3194c.gbf'
#genome_file = 'short_sequence.gb'

def match(sequence, index):
	print('\nmatch')
	print(index, sequence)

def get_aminoAcid(nucleotide_sequence, amino_sequence, start_position, stop_position, translation_table, strand_orientation):
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

#ask user for type (CDS, gene, tRNA)
outName = genome_file.replace('.gb', '_analysisOutput.txt')

with open(outName, 'w') as outputFile:
	outputFile.write('\t'.join(['Accession', 'Locus_Tag', 'Gene', 'Identity', 'RFP', 'index']) + '\n')  #RFP - reading frame position
	for record in SeqIO.parse(genome_file, 'genbank'):
		count = 1
	#	print(record.id)
		ID = record.id
		sequence = record.seq
		for feature in record.features:
			if feature.type=='CDS':
	#			print(feature)
				count+=1
				feats = 'CDS'
				start = feature.location.start.position
				end = feature.location.end.position
				strand = feature.location.strand
				locus_tag = feature.qualifiers.get('locus_tag')
				frame = feature.qualifiers.get('codon_start')
				gene = feature.qualifiers.get('gene')
				if type(gene) == list:
					gene = ' '.join(feature.qualifiers.get('gene'))
				transl_table = int(''.join(feature.qualifiers.get('transl_table')))
				AA_sequence = Seq(''.join(feature.qualifiers.get('translation')))
				extracted_feature = feature.extract(sequence).translate(table=transl_table)
				identity = 0
				
				if extracted_feature == AA_sequence:
					identity = 100.00
#					string = '\t'.join([str(ID), str(''.join(locus_tag)), str(gene), str(identity), str(''.join(frame)), str(1)])
#					outputFile.write(string + '\n')

				elif extracted_feature[-1] == '*':
					extracted_feature = extracted_feature[0:-1]
					if extracted_feature == AA_sequence:
						identity = 100.00
					elif 'M' + extracted_feature[1:] == AA_sequence:
						extracted_feature = 'M' + extracted_feature[1:-1]
						identity = round(100*len(extracted_feature)/float(len(AA_sequence)), 2)
					else:
						identity = 'NA'
#THIS SHOULD GO TO ANOTHER PROGRAM (MAUVE??) TO GET IDENTITY

				string = '\t'.join([str(ID), str(''.join(locus_tag)), str(gene), str(identity), str(''.join(frame)), str(1)])
				outputFile.write(string + '\n')

		print(count)
#start  = 49163 
#stop = 50149 
#AA = 'MDNTQKYAAIDLKSFYASVECILRKLDPLNTNLVVADESRTEKTICLAVSPALRSYNISGRLRLFELIQKVKTINYERLKIAKYFSAKSYNHLELINNPNLELDYIVAKPRMSTYIDYSSKIYSIYLKYFDPKDIHIYSIDEVFIDLTPYIKHYKLSADKLIENILFEILKTTQITATAGIGTNLYLAKIAMDILAKKQNINKDGLCIGYLDEMLYRRKLWQHTPINDFWRIGKGYATKLKSIGINNMGDLARYSLNNEDKLYQIFGVNTELLIDHAWGFESCTMQAIKEYKSKHISKVMAKVLPKPYSFKKARNMLKEIVDHMVRAN'
#get_aminoAcid(sequence, AA, start, stop)
