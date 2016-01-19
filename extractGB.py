'''
Created:	E. Reichenberger
Date:	10.30.2011
Modification Date: 10.23.2015

Purpose: To extact information from GenBank file and place the information to a new csv file. Output file should have the following format:
'Filename', 'Strain', 'DNA_Source', 'Locus_Tag', 'Product', 'Protein_ID', 'Note', 'Transl_Tbl', 'Seq_AA' (sep. by tabs)
NOTE: Some files do not have a complete genbank file (no translations or ORIGIN sequence). This will need to be addressed with fasta files
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


###############################################File List & Searchable Qualifiers###################################################
fileList = ['Campy1147c/Campy1147c_Chrom.gbf', 'Campy1147q/Campy1147q_Chrom_1.gbf', 'Campy1147q/Campy1147q_Chrom_2.gbf', 'Campy1147q/Campy1147q_Chrom_3.gbf', 'Campy1188c/Campy1188c_Chrom.gbf', 'Campy1188c/Campy1188c_Plasmid.gbf', 'Campy1246c/Campy1246c_Chrom.gbf', 'Campy1246c/Campy1246c_Plasmid.gbf', 'Campy1285c/Campy1285c_Chrom.gbf', 'Campy14076c/Campy14076c_Chrom.gbf', 'Campy3194c/Campy3194c_Chrom.gbf', 'Campy3194c/Campy3194c_Plasmid.gbf', 'CampyCP006702/CP006702_Chrom.gb', 'CampyCP007179/CP007179_Chrom.gb', 'CampyCP007181/CP007181_Chrom.gb']
qualies = ['locus_tag', 'product', 'protein_id', 'note', 'transl_table', 'translation']
#---------------------------------------------------------------------------------------------------------------------------------


############################Searchable Genbank Files & write qualifying features to output file####################################
for f in fileList:
	extension = os.path.splitext(f)[1]
	nameOut = ''
	if extension == '.gbf':
		nameOut = f.replace('gbf', 'txt')
	else:
		nameOut = f.replace('gb', 'txt') #not all files are gbf extension. Be sure to keep after above as replacing gb <- txtf
	strain = f.split('/')[0]
	file_name = f.split('/')[1].replace('.gbf', '')
	file_name = f.split('/')[1].replace('.gb', '')
	dna_type = f.split('_')[1]
	print(file_name)
	
	if dna_type.startswith('C'):
		dna_type = 'Chromosome'
	else:
		dna_type = 'Plasmid'
	
	string = ''
	with open(f, 'r') as handle:
		with open(nameOut, 'w') as outputFile:
			outputFile.write('\t'.join(['Filename', 'Strain', 'DNA_Source', 'Locus_Tag', 'Product', 'Protein_ID', 'Note', 'Transl_Tbl', 'Seq_AA']) + '\n')

			for record in SeqIO.parse(handle, 'genbank'):
				print(len(record.seq))
				for feature in record.features:
					if feature.type=='CDS':
						cds_dic = {} 
						string = '\t'.join([file_name, strain, dna_type])

						if f not in cds_dic:
							cds_dic[f] = {'filename':f, 'strain':strain, 'DNA_Source':dna_type, 'locus_tag':'NA', 'product':'NA', 'protein_id':'NA', 'note':'NA', 'transl_table':'NA', 'translation':'NA'}
						for q in qualies:
							cds_dic[f][q] = str(feature.qualifiers.get(q))

						for q in qualies:
							string = string + '\t' + strip_it(cds_dic[f][q])
						outputFile.write(string + '\n') # print(feature, string)
#---------------------------------------------------------------------------------------------------------------------------------

#This creates the corpus file
############################Create Corpus File (from individual text files)####################################
with open('Output/campy6_corpus_cds.txt', 'w') as outputFile:
	outputFile.write('\t'.join(['Filename', 'Strain', 'DNA_Source', 'Locus_Tag', 'Product', 'Protein_ID', 'Note', 'Transl_Tbl', 'Seq_AA']) + '\n')
	fileList = ['Campy1147c/Campy1147c_Chrom.txt', 'Campy1147q/Campy1147q_Chrom_1.txt', 'Campy1147q/Campy1147q_Chrom_2.txt', 'Campy1147q/Campy1147q_Chrom_3.txt', 'Campy1188c/Campy1188c_Chrom.txt', 'Campy1188c/Campy1188c_Plasmid.txt', 'Campy1246c/Campy1246c_Chrom.txt', 'Campy1246c/Campy1246c_Plasmid.txt', 'Campy1285c/Campy1285c_Chrom.txt', 'Campy14076c/Campy14076c_Chrom.txt', 'Campy3194c/Campy3194c_Chrom.txt', 'Campy3194c/Campy3194c_Plasmid.txt', 'CampyCP006702/CP006702_Chrom.txt', 'CampyCP007179/CP007179_Chrom.txt', 'CampyCP007181/CP007181_Chrom.txt']
	for f in fileList:
		with open(f, 'r') as inputFile:
			inputFile.readline() #skip writing first line
			lines = inputFile.readlines()
			for line in lines:	
				outputFile.write(line)
#---------------------------------------------------------------------------------------------------------------------------------
