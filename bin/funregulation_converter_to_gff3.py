"""
###########################################################
#
# funregulation_converter_to_gff3.py
#
# FunRegulation: Converter To GFF3 (Alexandre Rafael Lenz)
# Universidade do Estado da Bahia (UNEB) / Universidade de Caxias do Sul (UCS)
# https://alexandrelenz@github.com/alexandrelenz/funregulation.git
# Last update: 06.05.2020 (Lenz, A. R.)
#
# Gene regulatory networks (GRN) of 
# Penicillium ucsensis 2HH and Penicillium oxalicum 114-2 
# inferred by a computional biology approach
#
# We propose the inference of global GRNs for Penicillium ucsensis
# 2HH and Penicillium oxalicum 114-2, based on TF-TG orthology relationships of
# related species combined with TFBSs prediction. First, global GRNs of related
# species (A. nidulans, N. crassa and S. cerevisiae) afford the mapping of
# orthologous interactions. Further, the TFBSs prediction provides accuracy to
# TF-TG relationships.
#
# Regulatory sequences (promoter regions) of P. ucsensis 2HH and P. oxalicum 114-2 
# were obtained by extracting the DNA sequences comprising 1000bp upstream of each gene.
# This script uses annotation in gff3 format and whole genome sequences.
#
# Promoter sequences extraction requires annotation in GFF3 format.
# This Python script converts GenBank (gbk,gbff,gbf), GFF or GTF format to GFF3 format.
#
# # Python packages required:
# i)   pip3 install Biopython
#
#   INPUT FILE:
#
#    i) in_file_annotation (GenBank (gbk,gbff,gbf), GFF or GTF format)
#    
#   OUTPUT FILE:
#    
#    i) out_file_annotation (gff3 format)
#
###########################################################
"""

from BCBio import GFF
from Bio import SeqIO

in_file_type = 'genbank' # options: 'genbank', 'gbk', 'gbff', 'gbf', 'gff', 'gtf'
in_file_annotation = "/Users/arlenz/NetBeansProjects/funregulation/input/GRN/genomes/GCF_000182925.2_NC12_genomic.gbff"
out_file_annotation = "/Users/arlenz/NetBeansProjects/funregulation/input/GRN/genomes/GCF_000182925.2_NC12_genomic.gff3"

"""
    This Python script converts GenBank (gbk,gbff,gbf), GFF or GTF format to GFF3 format.
""" 
def convert_to_gff3():
    if in_file_type == 'genbank' or in_file_type == 'gbff' or in_file_type == 'gbf' or in_file_type == 'gbk':
        #Genbank format
        with open(out_file_annotation, "w") as out_handle, open(in_file_annotation) as in_handle:
            GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
    elif in_file_type == 'gff' or in_file_type == 'gtf':
        #GFF or GTF format
        with open(out_file_annotation, "w") as out_handle, open(in_file_annotation) as in_handle:
            GFF.write(GFF.parse(in_handle), out_handle)

    in_handle.close()
    out_handle.close()

""" 
    #################### MAIN ######################
"""
"""
    Main function of this program
""" 
if __name__ == '__main__':
    convert_to_gff3()