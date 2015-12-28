'''
Created:	E. Reichenberger
Date:	11.20.2015

Purpose: This script will open/parse genbank files, extract CDS features, and compare (e.g., calculate sequence identity) the amino acid translation to the sliced & translated nucleotide sequence. If a match or close match, the calculated identity, locus_tag, and accession number are written to *_analysisOutput.txt. If a match is not close, the nucleotide slice is written to a fasta file (*_unpairedAnnotations.fasta) where the header is comprised of '>' + locus_tag.

Input: Genbank file(s)
Output: *_analysisOutput.txt, *_unpairedAnnotations.fasta

Working python version: Python 3.4.3 :: Anaconda 2.3.0 (64-bit) :: Biopython 1.6
'''

import Bio
import os
from Bio import GenBank
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq

##########################################################################
genome_file = 'Campy3194c/Campy3194c.gb'

#ask user for type (CDS, product, tRNA) in future
outName = genome_file.split('/')[1].replace('.gb', '_analysisOutput.txt')
unpaired = genome_file.split('/')[1].replace('.gb', '_unpairedAnnotations.fasta')
#--------------------------------------------------------------------------

with open(unpaired, 'w') as unpairedFile: #fasta file for alignment
	with open(outName, 'w') as outputFile: #output file is input for R file
		outputFile.write('\t'.join(['Accession', 'Locus_Tag', 'Gene_Product', 'Identity', 'RFP', 'index']) + '\n')  #RFP - reading frame position

		##########################GET CDS INFORMATION#################################
		for record in SeqIO.parse(genome_file, 'genbank'):
			count = 1 #keeping track of CDS encounters
			ID = record.id
			sequence = record.seq
			for feature in record.features:
				if feature.type=='CDS':
					count+=1
					feats = 'CDS'
					locus_tag = feature.qualifiers.get('locus_tag')
					frame = feature.qualifiers.get('codon_start')
					product = feature.qualifiers.get('product')
					if type(product) == list:
						product = ' '.join(feature.qualifiers.get('product'))
					transl_table = int(''.join(feature.qualifiers.get('transl_table')))
					AA_sequence = Seq(''.join(feature.qualifiers.get('translation')))
					extracted_feature = feature.extract(sequence).translate(table=transl_table)
					identity = 0
					
					if extracted_feature == AA_sequence:
						identity = 100.00

					elif extracted_feature[-1] == '*':
						extracted_feature = extracted_feature[0:-1]
						if extracted_feature == AA_sequence:
							identity = 100.00
						elif 'M' + extracted_feature[1:] == AA_sequence:
							extracted_feature = 'M' + extracted_feature[1:-1]
							identity = round(100*len(extracted_feature)/float(len(AA_sequence)), 2)
						else: #THIS SHOULD GO TO ANOTHER PROGRAM (MAUVE??) TO GET IDENTITY
							identity = 'NA'
							fasta_entry = '>' + locus_tag + '\n' + feature.extract(sequence) + '\n'
							unpairedFile.write(fasta_entry)

					string = '\t'.join([str(ID), str(''.join(locus_tag)), str(product), str(identity), str(''.join(frame)), str(1)])
					outputFile.write(string + '\n')

			print(count)
#--------------------------------------------------------------------------
