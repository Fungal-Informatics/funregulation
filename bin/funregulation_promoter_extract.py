"""
###########################################################
#
# funregulation_promoter_extract.py
#
# FunRegulation: Promoter Extract (Alexandre Rafael Lenz)
# Universidade do Estado da Bahia (UNEB) / Universidade de Caxias do Sul (UCS)
# https://alexandrelenz@github.com/alexandrelenz/funregulation.git
# Last update: 06.05.2020 (Lenz, A. R.)
#
# Gene regulatory networks (GRN) of 
# Penicillium ucsensis 2HH and Penicillium oxalicum 114-2 
# inferred by a computational biology approach
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
# Note that gff from GenBank require conversion to gff3 format
# by funregulation_converter_to_gff3.py
#
#   INPUT FILES:
#
#    # Penicillium ucsensis 2HH (this study) genomic data:
#    i)  in_file_genome (fasta format)
#    ii) in_file_annotation (gff3 format)
#    # Penicillium oxalicum 114-2 (GCA000346795.1 pdev1.0 from GenBank) genomic data:
#    i)  in_file_genome (fasta format)
#    ii) in_file_annotation (gff3 format)
#    
#   OUTPUT FILES:
#    
#    # Tab-delimited output file:
#    i) out/extract_results/
#    # Promoter sequences fasta files:
#    i) out/extract_misc/
#
###########################################################
"""

import sys
import os
import inspect
import shutil
import urllib.parse
import lib.library as lib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from collections import namedtuple

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

"""
#################### INPUT PATHS ######################
"""
in_file_genome = '/Users/arlenz/NetBeansProjects/funregulation/GRN/input/genomes/GCA_000346795.1_pde_v1.0_genomic.fna'
in_file_annotation = '/Users/arlenz/NetBeansProjects/funregulation/input/Genomes/GCA_000346795.1_pde_v1.0_genomic.gff3'
out_folder = '/Users/arlenz/NetBeansProjects/funregulation/GRN/input/promoters/Poxalicum/'
out_file_promoters = os.path.join(out_folder,'extract_results','promoters.tsv')

"""
################ PROMOTER LENGHT DEFINITION ##############
"""
upstream = -1000
downstream = 0

"""
##############   NAMEDTUPLE  DEFINITIONS  #################
"""
# Initialized GeneInfo named tuple to handle with GFF3 annotation. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "ltype", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)

# Initialized PromoterInfo named tuple to handle with Promoter BED file. Note: namedtuple is immutable
PromoterInfoFields = ["gene", "strand", "scaffold", "start", "stop"]
PromoterRecord = namedtuple("PromoterRecord", PromoterInfoFields)

"""
#################### PRE-PROCESSING ######################
"""
#create folder structure
if not os.path.isdir(out_folder):
    os.makedirs(out_folder)
    os.makedirs(os.path.join(out_folder, 'extract_misc'))
    os.makedirs(os.path.join(out_folder, 'extract_results'))
    os.makedirs(os.path.join(out_folder, 'logfiles'))
else:
    if os.path.isdir(os.path.join(out_folder, 'extract_results')):
        shutil.rmtree(os.path.join(out_folder, 'extract_results'))
        os.makedirs(os.path.join(out_folder, 'extract_results'))
    #make sure subdirectories exist
    dirs = [os.path.join(out_folder, 'extract_misc'), os.path.join(out_folder, 'extract_results'), os.path.join(out_folder, 'logfiles')]
    for d in dirs:
        if not os.path.isdir(d):
            os.makedirs(d)
        
#create log file
log_name = os.path.join(out_folder, 'logfiles', 'funregulation-promoterextract.log')
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
#lib.SystemInfo()

#check input files to make sure they are not empty
input_checks = [in_file_genome, in_file_annotation]
for i in input_checks:
    if i:
        lib.checkinputs(i)

"""
#################### PROCESSING ######################
"""

"""
    Parse the GFF3 attribute column and return a dict
"""
def parse_gff_attributes(attributeString):
    if attributeString == ".": return {}
    ret ={}
    if ";" not in attributeString: 
        ret=attributeString
    if ";" in attributeString:
        for attribute in attributeString.split(";"):
            key, value = attribute.split("=")
            ret[urllib.parse.unquote(key)] = urllib.parse.unquote(value)
    return ret

"""
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.
"""
def parse_gff3(filename):
    lib.log.info("Parsing GFF3 File")
    with open(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else urllib.parse.unquote(parts[0]),
                "source": None if parts[1] == "." else urllib.parse.unquote(parts[1]),
                "ltype": None if parts[2] == "." else urllib.parse.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.parse.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.parse.unquote(parts[7]),
                "attributes": parse_gff_attributes(parts[8])
            }
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield GFFRecord(**normalizedInfo)
    lib.log.info("GFF3 File parsed correctly")
    
"""
    Write promoters information data (gene, scaffold, strand and start codon postions)
    in a tab-delimited file
"""
def write_promoters_information(out_file_promoters):
    lib.log.info("Creating promoters info file")
    recordCount = 0
    promoters_partially_extracted = 0
    
    if os.path.isfile(out_file_promoters):
        lib.log.info("File will be overwritten")
    else:
        lib.log.info("File not found. New promoters info file will be created.")
        f = open(os.path.join(out_folder, out_file_promoters), 'w')
        f.write('')
        f.close()
    
    with open(out_file_promoters, 'w') as outfp:
        record_list = list()
        for record in parse_gff3(in_file_annotation):
            record_list.append(record)
        last_end_position = 0
        last_seqid = ''
        last_strand = ''
        next_start_position = 0
        next_seqid = ''
        next_strand = ''
        pos = 0
        while (pos<len(record_list)):
            record = record_list[pos]
            if record.ltype == 'source':
                source_size = record.end
            if record.ltype == 'mRNA':
                #Access attributes like this: my_strand = record.strand
                #gene = str(record.attributes)
                
                if (pos+2 < len(record_list)):
                    next_start_position = record_list[pos+2].start
                    next_seqid = record_list[pos+2].seqid
                    next_strand = record_list[pos+2].strand
                
                if record.strand == '+':
                    if last_seqid != record.seqid or last_strand == '-' :
                        last_end_position = 0
                    if record.start-last_end_position+upstream > 0 :
                        outfp.write(record.attributes.get("locus_tag")+'\t'+record.strand+'\t'+record.seqid+'\t'+str(record.start+upstream)+'\t'+str(record.start+downstream)+'\n')
                    else:
                        # incomplete promoters
                        outfp.write(record.attributes.get("locus_tag")+'\t'+record.strand+'\t'+record.seqid+'\t'+str(last_end_position+1)+'\t'+str(record.start+downstream)+'\n')
                        lib.log.info("Promoter of gene " + record.attributes.get("locus_tag") + " can't be fully indentified")
                        promoters_partially_extracted += 1
                else:
                    if next_seqid != record.seqid:
                           next_start_position = source_size
                    if next_strand == '+':
                           next_start_position = 0
                    if next_strand == '+' or record.end-upstream-1-next_start_position < 0 :
                            outfp.write(record.attributes.get("locus_tag")+'\t'+record.strand+'\t'+record.seqid+'\t'+str(record.end-upstream-1)+'\t'+str(record.end-downstream-1)+'\n')
                    else:
                        # incomplete promoters
                        outfp.write(record.attributes.get("locus_tag")+'\t'+record.strand+'\t'+record.seqid+'\t'+str(next_start_position-1)+'\t'+str(record.end-downstream-1)+'\n')
                        lib.log.info("Promoter of gene " + record.attributes.get("locus_tag") + " can't be fully indentified")
                        promoters_partially_extracted += 1
                last_strand = record.strand
                last_end_position = record.end
                last_seqid = record.seqid
                recordCount += 1
            pos=pos+1
    
    lib.log.info("%d promoters were found" % recordCount)    
    lib.log.info("Promoters partially identified: %d" % promoters_partially_extracted)
    lib.log.info("Promoters info file successfully created")
    
"""
    Parse promoter information data (gene, scaffold, strand and start codon postions)
    from a tab-delimited file
"""
def parse_promoter_info_file(out_file_promoters):
    lib.log.info("Parsing promoters info file")
    with open(out_file_promoters) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(PromoterInfoFields)
            #Normalize data
            normalizedInfo = {
                "gene": None if parts[0] == "." else urllib.parse.unquote(parts[0]),
                "strand": None if parts[1] == "." else urllib.parse.unquote(parts[1]),
                "scaffold": None if parts[2] == "." else urllib.parse.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "stop": None if parts[4] == "." else int(parts[4])
            }
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield PromoterRecord(**normalizedInfo)
    lib.log.info("Promoters info file parsed correctly")

"""
    Extract promoter region sequences to a new fasta file
"""
def extract_promoters(out_file_promoters):
    lib.log.info("Start promoters extraction")
    promoters_partially_extracted = 0
    # parse fasta file and turn into dictionary
    records = SeqIO.to_dict(SeqIO.parse(open(in_file_genome), 'fasta'))

    for gene in parse_promoter_info_file(out_file_promoters):
        long_seq_record = records[gene.scaffold]
        long_seq = long_seq_record.seq
        alphabet = long_seq.alphabet
        short_seq = str(long_seq)[gene.start:gene.stop]
        if gene.strand == '-':
            short_seq = str(long_seq)[gene.stop:gene.start]
            my_dna = Seq(short_seq, generic_dna)
            my_dna = my_dna.reverse_complement()
            short_seq=str(my_dna)
            if gene.stop > len(my_dna)+gene.start :
                lib.log.info("Promoter of gene " + gene.gene + " can't be extracted")
                promoters_partially_extracted += 1
                
        short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=gene.gene, description=gene.scaffold+'_'+str(gene.start)+':'+str(gene.stop))
        with open(out_folder+'extract_misc/'+gene.gene+'_promoter.fasta', 'w') as f:
            SeqIO.write(short_seq_record, f, 'fasta')
    lib.log.info("Promoters partially extracted: %d" % promoters_partially_extracted)
    lib.log.info("Promoters in"+out_folder+'extract_misc')
    lib.log.info("Promoters successfully extracted")
    sys.exit(0)
    
""" 
    #################### MAIN ######################
"""
"""
    Main function of this program
""" 
if __name__ == '__main__':
    #Execute promoters info file writer
    write_promoters_information (out_file_promoters)
    #Execute promoter sequences extraction
    extract_promoters(out_file_promoters)