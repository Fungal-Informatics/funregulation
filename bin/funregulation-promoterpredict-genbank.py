#!/usr/bin/env python

import sys
import os
import inspect
import shutil
import argparse
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


#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='funregulation-promoterpredict.py', usage="%(prog)s [options] -g genome.fasta -a annotation.gff3 -o /output_folder",
    description = '''Script that does it all.''',
    epilog = "Written by Alexandre Lenz (2017-2019) arlenz@ucs.br / alenz@uneb.br",
    formatter_class = MyFormatter)
parser.add_argument('-g', '--genome', required=False, help='Genome in FASTA format')
parser.add_argument('-a', '--annotation', required=False, help='Annotation in GFF3 format')
parser.add_argument('-o', '--out', required=False, help='Output folder')
parser.add_argument('-s', '--species', required=False, help='Species name (e.g. "Aspergillus oryzae") use quotes if there is a space')
parser.add_argument('--strain', default=False, help='Strain name (e.g. CEA10)')
parser.add_argument('--name', required=False, help='Gene Name attribute from GFF3 - Default:Name')
parser.add_argument('--upstream', required=False, help='Promoter Upstream Lenght to be extracted - Default:-300')
parser.add_argument('--downstream', required=False, help='Promoter Downstream Lenght to be extracted - Default:+40')
args = parser.parse_args()

#testing without args
args.genome = '/Users/arlenz/NetBeansProjects/funregulation/input/N.crassa-OR74A.fasta'
args.annotation = '/Users/arlenz/NetBeansProjects/funregulation/input/N.crassa-OR74A.gff'
args.out = '/Users/arlenz/NetBeansProjects/funregulation/output/'
args.species = "Neurospora crassa"
args.strain = 'OR74A'
args.name = 'locus_tag'
args.upstream = -300
args.downstream = +40

#create folder structure
args.out = os.path.join(args.out,args.species,args.strain)
if not os.path.isdir(args.out):
    os.makedirs(args.out)
    os.makedirs(os.path.join(args.out, 'predict_misc'))
    os.makedirs(os.path.join(args.out, 'predict_results'))
    os.makedirs(os.path.join(args.out, 'logfiles'))
else:
    if os.path.isdir(os.path.join(args.out, 'predict_results')):
        shutil.rmtree(os.path.join(args.out, 'predict_results'))
        os.makedirs(os.path.join(args.out, 'predict_results'))
    #make sure subdirectories exist
    dirs = [os.path.join(args.out, 'predict_misc'), os.path.join(args.out, 'predict_results'), os.path.join(args.out, 'logfiles')]
    for d in dirs:
        if not os.path.isdir(d):
            os.makedirs(d)
        
#create log file
log_name = os.path.join(args.out, 'logfiles', 'funregulation-promoterpredict.log')
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
lib.SystemInfo()

#get version of funregulation
#version = lib.get_version()
#lib.log.info("Running %s" % version)

#check input files to make sure they are not empty
input_checks = [args.genome, args.annotation]
for i in input_checks:
    if i:
        lib.checkinputs(i)

# Initialized GeneInfo named tuple to handle with GFF3 annotation. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "ltype", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)

# Initialized PromoterInfo named tuple to handle with Promoter BED file. Note: namedtuple is immutable
bedInfoFields = ["gene", "strand", "scaffold", "start", "stop"]
BEDRecord = namedtuple("BEDRecord", bedInfoFields)

"""Parse the GFF3 attribute column and return a dict"""
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
Write promoters data (gene, scaffold, strand and start codon postions) to new bed file
 
"""
def create_promoters_bed_file(file_promoters_bed):
    lib.log.info("Start promoter prediction")
    recordCount = 0
    promoters_partially_extracted = 0
    
    if os.path.isfile(file_promoters_bed):
        lib.log.info("File does exist and will be overwritten")
    else:
        lib.log.info("No such file. New promoter prediction file will be created")
        f = open(os.path.join(args.out, file_promoters_bed), 'w')
        f.write('')
        f.close()
        
    
    with open(file_promoters_bed, 'w') as outfp:
        for record in parse_gff3(args.annotation):
            if record.ltype == 'gene':
                #Access attributes like this: my_strand = record.strand
                attributes = record.attributes
                gene = attributes[args.name]
                if record.strand == '+':
                    if record.start+args.upstream > 0:
                        outfp.write(gene+'\t'+record.strand+'\t'+record.seqid+'\t'+str(record.start+args.upstream)+'\t'+str(record.start+args.downstream+2)+'\n')
                    else:    
                        # incomplete promoters
                        outfp.write(gene+'\t'+record.strand+'\t'+record.seqid+'\t'+str(1)+'\t'+str(record.start+args.downstream+2)+'\n')
                        lib.log.info("Promoter from gene " + gene + " can't be fully extracted")
                        promoters_partially_extracted += 1
                else:
                    if record.end+args.upstream > 0:
                        outfp.write(gene+'\t'+record.strand+'\t'+record.seqid+'\t'+str(record.end+args.upstream)+'\t'+str(record.end+args.downstream-2)+'\n')
                    else:    
                        # incomplete promoters
                        outfp.write(gene+'\t'+record.strand+'\t'+record.seqid+'\t'+str(1)+'\t'+str(record.end+args.downstream-2)+'\n')
                        lib.log.info("Promoter from gene " + gene + " can't be fully extracted")
                        promoters_partially_extracted += 1
                recordCount += 1
    lib.log.info("%d promoters were found" % recordCount)
    lib.log.info("Promoters partially extracted: %d" % promoters_partially_extracted)
    lib.log.info("Promoter prediction successful")
    
"""
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.
    
"""
def parse_bed(file_promoters_bed):
    lib.log.info("Parsing BED File")

    with open(file_promoters_bed) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(bedInfoFields)
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
            yield BEDRecord(**normalizedInfo)
    lib.log.info("BED File parsed correctly")

"""
Reads promoters data (genes, scaffolds and postions) from bed file 
and extract each gene promoter to new fasta file
 
"""
def extract_promoters_from_genome_fasta_file(file_promoters_bed):
    lib.log.info("Start promoters extraction from genome")
    promoters_partially_extracted = 0
    # parse fasta file and turn into dictionary
    records = SeqIO.to_dict(SeqIO.parse(open(args.genome), 'fasta'))

    for gene in parse_bed(file_promoters_bed):
        long_seq_record = records[gene.scaffold]
        long_seq = long_seq_record.seq
        alphabet = long_seq.alphabet
        short_seq = str(long_seq)[gene.start:gene.stop]
        if gene.strand == '-':
            my_dna = Seq(short_seq, generic_dna)
            my_dna.reverse_complement()
            my_dna.complement()
            short_seq=str(my_dna)
            if gene.stop > len(my_dna)+gene.start :
                lib.log.info("Promoter from gene " + gene.gene + " can't be fully extracted")
                promoters_partially_extracted += 1
                
        short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=gene.gene, description=gene.scaffold+'_'+str(gene.start)+':'+str(gene.stop))
        with open(args.out+'/predict_misc/'+gene.gene+'_promoter.fasta', 'w') as f:
            SeqIO.write(short_seq_record, f, 'fasta')
    lib.log.info("Promoters partially extracted: %d" % promoters_partially_extracted)
    lib.log.info("Promoters extraction succesful")
    lib.log.info("Promoters available in "+args.out+'/predict_misc')
    lib.log.info("funregulation Promoter Prediction Successful")
    sys.exit(0)
    

# file definitions
file_promoters_bed = os.path.join(args.out,'predict_results','promoters.bed')
#Execute the promoters BED file creator
create_promoters_bed_file (file_promoters_bed)
#Execute the promoter extraction from Fasta genome
extract_promoters_from_genome_fasta_file(file_promoters_bed)