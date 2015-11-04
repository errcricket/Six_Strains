'''
Author:	E. Reichenberger
Date:		11.4.2015

Purpose:	Given a genome/plasmid nucleotide sequence, a gene as a protein sequence, & location of gene, 
extract the gene from the genome/plasmid nucleotide sequence.

CDS complement(49163..50149)
/note="DNA polymerase V subunit UmuC PRK03609; UV damage repair protein, ImpB/MucB/SamB family of Bacteria UniRef RepID=B9KDE6_CAMLR"
/transl_table=11
/product="UV damage repair protein, ImpB/MucB/SamB"


#to get compliment seq.complement(), to convert to AA seq.translate(table=2)
'''

import Bio
from Bio import SeqIO
from Bio import SeqFeature
from Bio.GenBank import Record
from Bio.GenBank.Record import Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqFeature import Reference, SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC

f_name = 'Campy3194c/Campy3194c_Plasmid.gbf'
record = SeqIO.parse(f_name, 'genbank').next() 
plasmid_sequence = record.seq

uv_gene_Psequence = 'MDNTQKYAAIDLKSFYASVECILRKLDPLNTNLVVADESRTEKTICLAVSPALRSYNISGRLRLFELIQKVKTINYERLKIAKYFSAKSYNHLELINNPNLELDYIVAKPRMSTYIDYSSKIYSIYLKYFDPKDIHIYSIDEVFIDLTPYIKHYKLSADKLIENILFEILKTTQITATAGIGTNLYLAKIAMDILAKKQNINKDGLCIGYLDEMLYRRKLWQHTPINDFWRIGKGYATKLKSIGINNMGDLARYSLNNEDKLYQIFGVNTELLIDHAWGFESCTMQAIKEYKSKHISKVMAKVLPKPYSFKKARNMLKEIVDHMVRAN'

start = 49163
stop = 50149

#uv_gene_Nsequence = plasmid_sequence[start-1:stop-1] #this is an off by one error, index does not start at zero??

uv_gene_Nsequence = plasmid_sequence[start:stop]
#complementN = uv_gene_Nsequence.complement()
RcomplementN = uv_gene_Nsequence.reverse_complement()

#TcomplementN = complementN.translate(table=11)
TRcomplementN = RcomplementN.translate(table=11)

#uv_gene_AAsequence = uv_gene_Nsequence.complement()
#uv_gene_AAsequence = uv_gene_Nsequence.translate(table=11)

#print(complementN[1:10])
#print(RcomplementN[1:10])
#print(TcomplementN[1:10])

if uv_gene_Psequence == TRcomplementN:
	print('match')

#print(uv_gene_Nsequence) #redirect output to fasta file.
