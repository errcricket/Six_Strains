'''
Created:	E. Reichenberger
Date:	10.30.2011
Modification Date: 10.23.2015

Purpose: To extact information from GenBank file and place the information to a new csv file. Output file should have the following format:
Strain\tBacteria\tStrain\tLocus_tag\tCDS_Region\tGene_info\tSequence\tProduct\n.
'''

import Bio
import os
import sys
from Bio import GenBank
from Bio import SeqIO
from Bio import SeqFeature
from Bio.GenBank import Record
from Bio.GenBank.Record import Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqFeature import Reference, SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
#from StringIO import StringIO

fileList = ['Campy1147c/Campy1147c_Chrom.gbf', 'Campy1147q/Campy1147q_Chrom_1.gbf', 'Campy1147q/Campy1147q_Chrom_2.gbf', 'Campy1147q/Campy1147q_Chrom_3.gbf', 'Campy1188c/Campy1188c_Chrom.gbf', 'Campy1188c/Campy1188c_Plasmid.gbf', 'Campy1246c/Campy1246c_Chrom.gbf', 'Campy1246c/Campy1246c_Plasmid.gbf', 'Campy1285c/Campy1285c_Chrom.gbf', 'Campy14076c/Campy14076c_Chrom.gbf', 'Campy3194c/Campy3194c_Chrom.gbf', 'Campy3194c/Campy3194c_Plasmid.gbf']
#fileList = ['Campy1147c/Campy1147c_Chrom.gb', 'Campy1147q/Campy1147q_Chrom_1.gb', 'Campy1147q/Campy1147q_Chrom_2.gb', 'Campy1147q/Campy1147q_Chrom_3.gb', 'Campy1188c/Campy1188c_Chrom.gb', 'Campy1188c/Campy1188c_Plasmid.gb', 'Campy1246c/Campy1246c_Chrom.gb', 'Campy1246c/Campy1246c_Plasmid.gb', 'Campy1285c/Campy1285c_Chrom.gb', 'Campy14076c/Campy14076c_Chrom.gb', 'Campy3194c/Campy3194c_Chrom.gb', 'Campy3194c/Campy3194c_Plasmid.gb']
#fileList = ['Campy1147c/Campy1147c.gb', 'Campy14076c/Campy14076c.gb', 'Campy3194c/Campy3194c.gb', 'Campy1246c/Campy1246c.gb', 'Campy1147q/Campy1147q.gb', 'Campy1285c/Campy1285c.gb', 'Campy1188c/Campy1188c.gb']
campy_dic = {}

fileCounter = 0 #will increment +1 for each loop & will be the unique number
uniqueID = str(fileCounter) #each record must have a unique ID
qualies = {'locus_tag':'NA', 'product':'NA', 'protein_id':'NA', 'note':'NA', 'transl_table':'NA', 'translation':'NA'}

def strip_it(string_name):
	stripers = ['[', ']', '\'', '"']
	for s in stripers:
		string_name = string_name.replace(s, '')
		string_name = string_name.lstrip()
	return string_name

for f in fileList:
	nameOut = f.replace('gbf', 'txt')
	strain = f.split('/')[0]
	file_name = f.split('/')[1].replace('.gbf', '')
	dna_type = f.split('_')[1]
	#print(file_name, dna_type)
	print(file_name)
	
	if dna_type.startswith('C'):
		dna_type = 'Chromosome'
	else:
		dna_type = 'Plasmid'
	
	with open(f, 'r') as inputFile:
		with open(nameOut, 'w') as outputFile:
			outputFile.write('Filename\tStrain\tDNA_Source\tLocus_Tag\tProduct\tTransl_Tbl\tNote\tSeq_AA\tProtein_ID\n')
			record = SeqIO.parse(f, 'genbank').next()

			featureCount = 0
			for f in record.features:
				if f.type == 'CDS':
					CDS = record.features[featureCount]
					string = file_name + '\t' + strain + '\t' + dna_type 

					for q in qualies:
						if q in CDS.qualifiers:
							string = string + '\t' + str(CDS.qualifiers[q]).replace('\'', '')
							qualies[q] = CDS.qualifiers[q]
						else: # not q in CDS.qualifiers:
							string = string + '\tNA' 

					string = strip_it(string)
					outputFile.write(string+'\n')
						
				featureCount+=1

#This creates the corpus file
with open('Output/campy6_corpus_cds.txt', 'w') as outputFile:
	outputFile.write('Filename\tStrain\tDNA_Source\tLocus_Tag\tProduct\tTransl_Tbl\tNote\tSeq_AA\tProtein_ID\n')
	fileList = ['Campy1147c/Campy1147c_Chrom.txt', 'Campy1147q/Campy1147q_Chrom_1.txt', 'Campy1147q/Campy1147q_Chrom_2.txt', 'Campy1147q/Campy1147q_Chrom_3.txt', 'Campy1188c/Campy1188c_Chrom.txt', 'Campy1188c/Campy1188c_Plasmid.txt', 'Campy1246c/Campy1246c_Chrom.txt', 'Campy1246c/Campy1246c_Plasmid.txt', 'Campy1285c/Campy1285c_Chrom.txt', 'Campy14076c/Campy14076c_Chrom.txt', 'Campy3194c/Campy3194c_Chrom.txt', 'Campy3194c/Campy3194c_Plasmid.txt']
	for f in fileList:
		with open(f, 'r') as inputFile:
			inputFile.readline()
			lines = inputFile.readlines()
			for line in lines:	
				outputFile.write(line)

