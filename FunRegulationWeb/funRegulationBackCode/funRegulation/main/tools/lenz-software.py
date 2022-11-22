"""
###########################################################
#
# load_data.py
#
# FunRegulation: Network Inference (Alexandre Rafael Lenz)
# Universidade do Estado da Bahia (UNEB)
# https://github.com/Fungal-Informatics/funregulation
# Last update: 21.09.2021 (Lenz, A. R.)
#
#   INPUT DATA:
#
# a) Genomic data:
#    Data were collected and organized in tab-delimited files.
#    i) Organisms
#    ii) Genes
#    iii) Proteins
#    iv) Annotations
#
# b) Cis-BP:
#    Data were collected and organized in tab-delimited files for each species.
#    i) TF and PWM data
#
###########################################################
"""
GRNInferenceVersion = '1.0 - 04/10/2021'

import psycopg2
import sys, re
import os, platform
import urllib.parse
import lib.library as lib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import namedtuple
from suds.client import Client

#Initialize Folder Paths
in_folder = "/Users/arlenz/Documents/Alexandre/Trabalho/Pesquisa/FunWorld/Tools/FunRegulation/Software/Database/"
out_folder = "/Users/arlenz/Documents/Alexandre/Trabalho/Pesquisa/FunWorld/Tools/FunRegulation/Software/Database/Results/"

#Initialize file Paths
in_file_organisms = os.path.join(in_folder,'Ensembl_Species.tsv')
in_file_genes = os.path.join(in_folder,'Model/Aspergillus_nidulans/A_nidulans_FGSC_A4_current_features.gff3')
in_file_proteins = os.path.join(in_folder,'Model/Aspergillus_nidulans/A_nidulans_FGSC_A4_current_orf_trans_all.fasta')
in_file_pwms = os.path.join(in_folder,'Cis-BP/Species/Model/Aspergillus_nidulans_2021_09_17_6_25_pm/TF_Information.txt')
in_file_model_regulatory = os.path.join(in_folder,'Model/Aspergillus_nidulans.tsv')
in_file_orthology = os.path.join(in_folder,'Model/Fusarium_graminearum/Anidulans.proteinortho.tsv')
in_file_genome = os.path.join(in_folder,'Model/Fusarium_graminearum/Fusarium_graminearum_gca_000240135.ASM24013v3.dna.toplevel.fa')
#Initialize Output File Paths
#...

#DB definitions
dbConnection = None;

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

#create log file
log_name = os.path.join(out_folder, 'logfiles', 'funregulation.log')
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

"""
    Class defining an Organism
"""
class Organism:
    def __init__(self, order, genus, species, strain, taxon, assembly_name, assembly_accession, is_model, cis_bp):
        self.order = order
        self.genus = genus
        self.species = species
        self.strain = strain
        self.taxon = taxon
        self.assembly_name = assembly_name
        self.assembly_accession = assembly_accession
        self.is_model = is_model
        self.cis_bp = cis_bp

"""
    Class defining a Gene
"""
class Gene:
    def __init__(self, organism, locus_tag, symbol_gene, description, is_tf):
        self.organism = organism
        self.locus_tag = locus_tag
        self.symbol_gene = symbol_gene
        self.description = description
        self.is_tf = is_tf
        
"""
    Class defining a Promoter
"""
class Promoter:
    def __init__(self, locus_tag, strand, source, start, stop):
        self.locus_tag = locus_tag
        self.strand = strand
        self.source = source
        self.start = start
        self.stop = stop

"""
    Class defining a Promoter
"""
class Protein:
    def __init__(self, locus_tag, id, interpro, pfam, go, gene3d, reactome, panther, uniprot, kegg_enzyme, cazy, uniparc):
        self.locus_tag = locus_tag
        self.id = id
        self.interpro = interpro
        self.pfam = pfam
        self.go = go
        self.gene3d = gene3d
        self.reactome = reactome
        self.panther = panther
        self.uniprot = uniprot
        self.kegg_enzyme = kegg_enzyme
        self.cazy = cazy
        self.uniparc = uniparc

"""
    Class defining an Orthology
"""
class Orthology:
    def __init__(self, model_protein, target_protein):
        self.model_protein = model_protein
        self.target_protein = target_protein
        
        
"""
    Class defining a PWM
"""
class Pwm:
    def __init__(self, id, locus_tag, motif_id, status, tf_family, motif_type, msource_author, msource, pubmedid):
        self.id = id
        self.locus_tag = locus_tag
        self.motif_id = motif_id
        self.status = status
        self.tf_family = tf_family
        self.motif_type = motif_type
        self.msource_author = msource_author
        self.msource = msource
        self.pubmedid = pubmedid
        
"""
    Class defining a Model Regulatory Interaction
"""
class ModelRegulatory:
    def __init__(self, id, tf_locus_tag, tg_locus_tag, regulatory_function, evidence, experiment, experimental_condition, pubmedid, publication):
        self.id = id
        self.tf_locus_tag = tf_locus_tag
        self.tg_locus_tag = tg_locus_tag
        self.regulatory_function = regulatory_function
        self.evidence = evidence
        self.experiment = experiment
        self.experimental_condition = experimental_condition
        self.pubmedid = pubmedid
        self.publication = publication

"""
    Class defining a new Regulatory Interaction
"""
class RegulatoryInteraction:
    def __init__(self, id, tf_locus_tag, tg_locus_tag, regulatory_function, pubmedid_source):
        self.id = id
        self.tf_locus_tag = tf_locus_tag
        self.tg_locus_tag = tg_locus_tag
        self.regulatory_function = regulatory_function
        self.pubmedid_source = pubmedid_source

"""
    Class defining a new TFBS prediction
"""
class TFBS:
    def __init__(self, id, regulatory_interaction_id, pwm_id, strand, start, end, sequence, weight, pval, ln_pval, sig):
        self.id = id
        self.regulatory_interaction_id = regulatory_interaction_id
        self.pwm_id = pwm_id
        self.strand = strand
        self.start = start
        self.end = end
        self.sequence = sequence
        self.weight = weight
        self.pval = pval
        self.ln_pval = ln_pval
        self.sig = sig

""" create a database connection to a SQLite database
"""
def create_db_connection():
    try:
        con = psycopg2.connect(host='localhost', database='funregulation',
        user='postgres', password='123456')
        lib.log.info("Successfully Connected to PostgreSQL")
        return con
    except (Exception, psycopg2.Error) as error:
        lib.log.info(error)
        

"""
    Insert organism
"""
def insert_organism(organism):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO organism VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                            (organism.order, 
                                                            organism.genus, 
                                                            organism.species, 
                                                            organism.strain, 
                                                            organism.taxon, 
                                                            organism.assembly_name, 
                                                            organism.assembly_accession, 
                                                            organism.is_model, 
                                                            organism.cis_bp))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE organism")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE organism ")
        lib.log.info(str(organism.order) + " " +
                    str(organism.genus) + " " + 
                    str(organism.species) + " " +
                    str(organism.strain))


"""
    Parse Ensembl Species input file
"""
def parse_organism_file(filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            order = urllib.parse.unquote(parts[0])
            genus = urllib.parse.unquote(parts[1])
            species = urllib.parse.unquote(parts[2])
            strain = urllib.parse.unquote(parts[3])
            taxon = urllib.parse.unquote(parts[4])
            assembly_name = urllib.parse.unquote(parts[5])
            assembly_accession = urllib.parse.unquote(parts[6])
            is_model = urllib.parse.unquote(parts[7])
            cis_bp = urllib.parse.unquote(parts[8])
            organism = Organism(order,genus,species,strain,taxon,assembly_name,assembly_accession,is_model,cis_bp)
            insert_organism(organism)
    in_file.close()
    lib.log.info(filename + " parsed correctly")

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
def parse_gff3_file(filename):
    lib.log.info("Parsing "+ filename)
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
    Select organism by assembly
"""
def select_organism_by_assembly_name(source):
    organism = 0
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT id FROM organism WHERE assembly_name = %s"
        cursor.execute(postgreSQL_select_Query, (source,))
        records = cursor.fetchall()
        for row in records:
            organism = row[0]
        return organism
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table organism", error)
        lib.log.info(source)

"""
    Insert gene
"""
def insert_gene(gene):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO gene VALUES (%s, %s, %s, %s, %s)",
                                                        (gene.organism, 
                                                        gene.locus_tag, 
                                                        gene.symbol_gene, 
                                                        gene.description, 
                                                        gene.is_tf))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE gene")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE gene", error)
        lib.log.info(str(gene.organism) + " " +
                     str(gene.locus_tag) + " " + 
                     str(gene.symbol_gene) + " " +
                     str(gene.description) + " " +
                     str(gene.is_tf))
                     
"""
    Update gene
"""
def update_gene(gene):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("UPDATE gene SET symbol_gene = %s, description = %s, is_tf = %s WHERE organism = %s AND locus_tag = %s",
                                (gene.symbol_gene, 
                                gene.description,
                                gene.is_tf,
                                gene.organism,
                                gene.locus_tag))
        dbConnection.commit()
        lib.log.info("Record updated successfully into TABLE gene")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to update data into TABLE gene", error)
        lib.log.info(str(gene.organism) + " " +
                     str(gene.locus_tag) + " " + 
                     str(gene.symbol_gene) + " " +
                     str(gene.description) + " " +
                     str(gene.is_tf))
                    
"""
    Insert promoter
"""
def insert_promoter(promoter):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO promoter VALUES (%s, %s, %s, %s, %s)",
                                                            (promoter.locus_tag, 
                                                            promoter.strand, 
                                                            promoter.source, 
                                                            promoter.start, 
                                                            promoter.stop))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE promoter")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE promoter", error)
        lib.log.info(str(promoter.locus_tag) + " " +
                    str(promoter.strand) + " " + 
                    str(promoter.source) + " " +
                    str(promoter.start) + " " +
                    str(promoter.stop))

"""
    GFF3 handler
"""
def gff3_handler(in_file_genes):
    lib.log.info("Parsing "+ in_file_genes)
    recordCount = 0
    promoters_partially_extracted = 0
    organism_id = 0
    
    record_list = list()
    for record in parse_gff3_file(in_file_genes):
        record_list.append(record)
        
    pos = 0
    while (pos<len(record_list)):
        record = record_list[pos]
        if record.ltype == 'chromosome' or record.ltype == 'supercontig':
            source_size = record.end
            source = record.source
            organism_id = select_organism_by_assembly_name(source)
        else:
            if (record.ltype == 'gene' or 
                record.ltype == 'pseudogene' or 
                record.ltype == 'transposable_element_gene' or 
                record.ltype == 'blocked_reading_frame'):
                #Access attributes like this: my_strand = record.strand
                #gene = str(record.attributes)
                locus_tag = record.attributes.get("ID")
                #description = record.attributes.get("description")
                description = None #record.attributes.get("description")
                symbol_gene = ''
                if record.attributes.get("Name") is not None:
                    symbol_gene = record.attributes.get("Name")
                is_tf = False
                gene = Gene(organism_id,locus_tag,symbol_gene,description,is_tf)
                insert_gene(gene)

                promoter = None
                if record.strand == '+':
                    if record.start+upstream > 0 :
                        promoter = Promoter(locus_tag, record.strand, record.seqid, record.start+upstream, record.start+downstream)
                    else:
                        # incomplete promoters
                        promoter = Promoter(locus_tag, record.strand, record.seqid, 1, record.start+downstream)
                        lib.log.info("Promoter of gene " + locus_tag + " can't be fully indentified")
                        promoters_partially_extracted += 1
                else:
                    if record.end-source_size <= 0 :
                        promoter = Promoter(locus_tag, record.strand, record.seqid, record.end-upstream, record.end-downstream)
                    else:
                        # incomplete promoters
                        promoter = Promoter(locus_tag, record.strand, record.seqid, source_size, record.end-downstream)
                        lib.log.info("Promoter of gene " + locus_tag + " can't be fully indentified")
                        promoters_partially_extracted += 1
                recordCount += 1
                insert_promoter(promoter)
        pos=pos+1
    
    lib.log.info("%d genes were found" % recordCount)
    lib.log.info("Promoters partially identified: %d" % promoters_partially_extracted)
    lib.log.info("GFF3 file successfully parsed")

"""
    Insert protein
"""
def insert_protein(protein):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO protein VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                            (protein.locus_tag, 
                                                            protein.id, 
                                                            protein.interpro, 
                                                            protein.pfam, 
                                                            protein.go,
                                                            protein.gene3d,
                                                            protein.reactome,
                                                            protein.panther,
                                                            protein.uniprot,
                                                            protein.kegg_enzyme,
                                                            protein.cazy,
                                                            protein.uniparc
                                                            ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE protein")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE protein", error)
        lib.log.info(str(protein.locus_tag) + " " +
                    str(protein.id))

"""
    Parse Protein Fasta file
"""
def parse_protein_file(in_file_proteins):
    lib.log.info("Parsing "+ in_file_proteins)
    for rec in SeqIO.parse(in_file_proteins, 'fasta'):
        
        #when locus_tag != protein_id
        #rec.description = re.search(r'gene:(.*?) transcript:', rec.description).group(1)
        #protein = Protein(rec.description, rec.id,'','','','','','','','','','')
        
        #when locus_tag == protein_id
        protein = Protein(rec.id, rec.id,'','','','','','','','','','')
        insert_protein(protein)
        
"""
    Insert pwm
"""
def insert_pwm(pwm):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO pwm VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                        (pwm.locus_tag, 
                                                        pwm.motif_id, 
                                                        pwm.status, 
                                                        pwm.tf_family, 
                                                        pwm.motif_type,
                                                        pwm.msource_author,
                                                        pwm.msource,
                                                        pwm.pubmedid
                                                        ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE pwm")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE pwm", error)
        lib.log.info(str(pwm.locus_tag) + " " +
                    str(pwm.motif_id) + " " +
                    str(pwm.status) + " " +
                    str(pwm.tf_family) + " " +
                    str(pwm.motif_type) + " " +
                    str(pwm.msource_author) + " " +
                    str(pwm.msource) + " " +
                    str(pwm.pubmedid))
           
"""
    Parse CIS-BP PWM input file
"""
def parse_pwm_file(in_file_pwm):
    lib.log.info("Parsing "+ in_file_pwm)
    with open(in_file_pwm) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            line_parts = line.strip().split("\t")
            
            motif_id = urllib.parse.unquote(line_parts[3]) 
            if motif_id != '.' and motif_id != 'Motif_ID':
                locus_tag = urllib.parse.unquote(line_parts[5])
                status = urllib.parse.unquote(line_parts[8])
                tf_family = urllib.parse.unquote(line_parts[9])
                motif_type = urllib.parse.unquote(line_parts[14])
                msource = urllib.parse.unquote(line_parts[16])
                msource_author = urllib.parse.unquote(line_parts[17])
                pubmedid = urllib.parse.unquote(line_parts[19])
                if pubmedid == 'NULL':
                    pubmedid = ''
                pwm = Pwm(0,locus_tag, motif_id, status, tf_family, motif_type, msource_author, msource, pubmedid)
                insert_pwm(pwm)
    in_file.close()
    lib.log.info(in_file_pwm + " parsed correctly")


"""
    Select gene by locus_tag
"""
def select_gene_by_locus_tag(locus_tag):
    gene = None
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * FROM gene WHERE locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (locus_tag,))
        rec = cursor.fetchall()
        for row in rec:
            gene = Gene(row[0],row[1],row[2],row[3],row[4])
            return gene
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table gene", error)

"""
    Insert Model regulatory interaction
"""
def insert_model_regulatory(model_regulatory):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO model_regulatory VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                                    (model_regulatory.tf_locus_tag, 
                                                                    model_regulatory.tg_locus_tag, 
                                                                    model_regulatory.regulatory_function, 
                                                                    model_regulatory.evidence, 
                                                                    model_regulatory.experiment,
                                                                    model_regulatory.experimental_condition,
                                                                    model_regulatory.pubmedid,
                                                                    model_regulatory.publication
                                                                    ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE model_regulatory")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE model_regulatory", error)
        lib.log.info(str(model_regulatory.tf_locus_tag) + " " +
                    str(model_regulatory.tg_locus_tag) + " " +
                    str(model_regulatory.regulatory_function))
                    
"""
    Parse Model regulatory interactions file
"""
def parse_model_regulatory_file(filename):
    lib.log.info("Parsing " + filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            
            #update gene as TF
            tf_locus_tag = urllib.parse.unquote(parts[0])
            tf = select_gene_by_locus_tag(tf_locus_tag)
            tf.is_tf = 'True'
            if tf.symbol_gene is None or tf.symbol_gene == '':
                tf.symbol_gene = urllib.parse.unquote(parts[1])
            update_gene(tf)
            
            #update gene as TG
            tg_locus_tag = urllib.parse.unquote(parts[2])
            tg = select_gene_by_locus_tag(tg_locus_tag)
            if tg is not None:
                if tg.symbol_gene is None or tg.symbol_gene == '':
                    tg.symbol_gene = urllib.parse.unquote(parts[3])
                update_gene(tg)
                regulatory_function = urllib.parse.unquote(parts[4])
                evidence = urllib.parse.unquote(parts[5])
                experiment = urllib.parse.unquote(parts[6])
                experimental_condition = urllib.parse.unquote(parts[7])
                pubmedid = urllib.parse.unquote(parts[8])
                publication = urllib.parse.unquote(parts[9])
                model_regulatory = ModelRegulatory(0,tf_locus_tag,tg_locus_tag,regulatory_function,evidence,experiment,experimental_condition,pubmedid,publication)
                insert_model_regulatory(model_regulatory)
    in_file.close()
    lib.log.info(filename + " parsed correctly")

"""
    Select protein by id
"""
def select_protein_by_id(protein_id):
    protein = None
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * FROM protein WHERE id = %s"
        cursor.execute(postgreSQL_select_Query, (protein_id,))
        rec = cursor.fetchall()
        for row in rec:
            protein = Protein(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11])
            return protein
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table protein", error)
        lib.log.info(source)

"""
    Insert Orthology
"""
def insert_orthology(orthology):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO orthology VALUES (%s, %s)",
                                                                    (orthology.model_protein.id, 
                                                                    orthology.target_protein.id))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE orthology")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE orthology", error)
        lib.log.info(str(orthology.model_protein.id) + " " +
                    str(orthology.target_protein.id))

"""
    Parse orthology file - ProteinOrtho
"""
def parse_orthology_file(filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            line_parts = line.strip().split("\t")
            
            model = urllib.parse.unquote(line_parts[3])
            target = urllib.parse.unquote(line_parts[4])
            
            model_parts = model.strip().split(",")
            target_parts = target.strip().split(",")
            
            for record_model in model_parts:
                for record_target in target_parts:
                    if (record_model != '*' and record_target != '*'):
                        model_protein = select_protein_by_id(record_model)
                        target_protein = select_protein_by_id(record_target)
                        orthology = Orthology(model_protein,target_protein)
                        insert_orthology(orthology)
                    
    in_file.close()
    lib.log.info(filename + " parsed correctly")

"""
    Select model regulatory interactions by organism id
"""
def select_model_regulatory_by_organism_id(organism_id):
    model_regulatory_interactions = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT DISTINCT * from model_regulatory model right join gene gen on model.tf_locus_tag = gen.locus_tag AND gen.organism = %s WHERE model.tf_locus_tag IS NOT NULL", (organism_id))
        rec = cursor.fetchall()
        for row in rec:
            model_regulatory = ModelRegulatory(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8])
            model_regulatory_interactions.append(model_regulatory)
        cursor.close()
        return model_regulatory_interactions
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table model_regulatory", error)
        lib.log.info(organism_id)

"""
    Select orthologs by model_locus_tag and target organism id
"""
def select_orthologs_by_target_organism(model_locus_tag, target_organism_id):
    orthology_list = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT DISTINCT model_protein, target_protein from orthology ortho join protein prot on ortho.model_protein = prot.id join gene gen on prot.locus_tag = %s AND gen.organism = %s WHERE ortho.model_protein IS NOT NULL", (model_locus_tag, target_organism_id))
        rec = cursor.fetchall()
        for row in rec:
            ortho = Orthology(select_protein_by_id(row[0]),select_protein_by_id(row[1]))
            orthology_list.append(ortho)
        cursor.close()
        return orthology_list
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table orthology", error)
        lib.log.info(model_locus_tag)

"""
    Select locus_tag by protein id
"""
def select_locus_tag_by_protein_id(protein_id):
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT locus_tag from protein WHERE id = %s"
        cursor.execute(postgreSQL_select_Query, (protein_id,))
        rec = cursor.fetchone()
        print(rec)
        cursor.close()
        return rec
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table protein", error)


"""
    Insert NEW regulatory interaction
"""
def insert_regulatory_interaction(regulatory_interaction):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO regulatory_interaction VALUES (default, %s, %s, %s, %s)",
                                                                    (regulatory_interaction.tf_locus_tag, 
                                                                    regulatory_interaction.tg_locus_tag, 
                                                                    regulatory_interaction.regulatory_function,
                                                                    regulatory_interaction.pubmedid_source
                                                                    ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE regulatory_interaction")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE regulatory_interaction", error)
        lib.log.info(str(regulatory_interaction.tf_locus_tag) + " " +
                    str(regulatory_interaction.tg_locus_tag) + " " +
                    str(regulatory_interaction.regulatory_function))
                    
"""
    Construct GRN Orthology
"""
def construct_grn_orthology(model_organism_id, target_organism_id):
    
    model_regulatory_interactions = select_model_regulatory_by_organism_id(model_organism_id)
    
    for model_regulatory in model_regulatory_interactions:
        tf_orthologs = select_orthologs_by_target_organism(model_regulatory.tf_locus_tag, target_organism_id)
        tg_orthologs = select_orthologs_by_target_organism(model_regulatory.tg_locus_tag, target_organism_id)
        if len(tf_orthologs)!=0 and len(tg_orthologs)!=0:
            for ortholog_tf in tf_orthologs:
                for ortholog_tg in tg_orthologs:
                    regulatory_interaction = RegulatoryInteraction(0,ortholog_tf.target_protein.locus_tag, ortholog_tg.target_protein.locus_tag, model_regulatory.regulatory_function,model_regulatory.pubmedid)
                    
                    insert_regulatory_interaction(regulatory_interaction)
                    
                    #update gene as TF
                    tf = select_gene_by_locus_tag(ortholog_tf.target_protein.locus_tag)
                    tf.is_tf = 'True'
                    update_gene(tf)


"""
##################### RSAT WS ########################
"""
""" 
    Create the SOAP client to request RSAT services
""" 
# Define URL for RSAT services (alternative RSAT WSs to use when the UNAM WS is down)
#wsdlUrl =  'http://rsat.ulb.ac.be/rsat/web_services/RSATWS.wsdl'
#wsdlUrl =  'http://rsat01.biologie.ens.fr/rsat/web_services/RSATWS.wsdl'
#wsdlUrl = 'http://pedagogix-tagc.univ-mrs.fr/rsat/web_services/RSATWS.wsdl'
wsdlUrl = 'http://embnet.ccg.unam.mx/rsat/web_services/RSATWS.wsdl'
#wsdlUrl = 'http://rsat-tagc.univ-mrs.fr/rsat/web_services/RSATWS.wsdl'

# Create the client
try:
    client = Client(wsdlUrl)
except (Exception) as error:
    lib.log.info("Connection to RSAT Server failed!", error)
    raise
# Need the service interface to perform requests
rsat_service = client.service

"""
    Define client header (optional)
"""
userAgent = 'RSAT-Client/v%s (%s; Python %s; %s)' % (
    GRNInferenceVersion, 
    os.path.basename( __file__ ),
    platform.python_version(), 
    platform.system()
)
httpHeaders = {'User-agent': userAgent}
client.set_options(headers=httpHeaders)
client.set_options(timeout=300)




"""
    RSAT Web Service: Scan sequences for a given matrix
"""
def call_matrix_scan(service, fasta_content_str, matrix_str):
    #RSAT config parameters
    uth_pval = '1e-2'
    matrix_format = "cis-bp"
    
    print(fasta_content_str)
    print(matrix_str)
    # Wrap all arguments into a named list
    #http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
    arguments = {
            'sequence' : fasta_content_str,
            'matrix' : matrix_str,
            'matrix_format' : matrix_format,
            'uth' : ['pval '+uth_pval],
            'lth' : ['score '+str(1)],
            #'quick' : 1,
            'str' : 2,
            'origin' : 'start',
            'background_input' : 1, # this option requires 'markov'
            'markov' : 1,
            'pseudo' : 1,
            'decimals' : 1,
            'n_treatment' : 'score',
            'background_pseudo': 0.01
            #'sequence_format' : 'fasta',
            #'organism': 'fungi'
            #'verbosity' : 1
    }

    ## Perform SOAP request on RSAT server
    result = service.matrix_scan(arguments)
    print (result.client)
    return result.client

"""
    Select tfs by organism id
"""
def select_tfs_by_organism(organism_id):
    tf_list = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT * from gene WHERE is_tf = 'True' and organism = %s", (organism_id))
        rec = cursor.fetchall()
        for row in rec:
            tf_list.append(Gene(row[0],row[1],row[2],row[3],row[4]))
        cursor.close()
        return tf_list
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table gene", error)
        lib.log.info(organism_id)

"""
    Select pwms by tf_locus_tag
"""
def select_pwms_by_locus_tag(locus_tag):
    pwm_list = list()
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * from pwm WHERE locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (locus_tag,))
        rec = cursor.fetchall()
        for row in rec:
            pwm_list.append(Pwm(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8]))
        cursor.close()
        return pwm_list
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table pwm", error)
        lib.log.info(locus_tag)

"""
    Select regulatory interactions by tf_locus_tag
"""
def select_regulatory_interactions_by_tf_locus_tag(tf_locus_tag):
    regulatory_interactions = list()
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * from regulatory_interaction WHERE tf_locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (tf_locus_tag,))
        rec = cursor.fetchall()
        for row in rec:
            regulatory_interactions.append(RegulatoryInteraction(row[0],row[1],row[2],row[3],row[4]))
        cursor.close()
        return regulatory_interactions
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table regulatory_interaction", error)
        lib.log.info(tf_locus_tag)

"""
    Select promoter by tg_locus_tag
"""
def select_promoter_by_locus_tag(tg_locus_tag):
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * from promoter WHERE locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (tg_locus_tag,))
        rec = cursor.fetchone()
        promoter = Promoter(rec[0],rec[1],rec[2],rec[3],rec[4])
        cursor.close()
        return promoter
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table promoter", error)
        lib.log.info(tg_locus_tag)

"""
    Extract promoter region sequence
"""
def extract_promoter(genome, tg_locus_tag):
    
    promoter = select_promoter_by_locus_tag(tg_locus_tag)
    
    long_seq_record = genome[promoter.source]
    long_seq = long_seq_record.seq
    short_seq = str(long_seq)[promoter.start:promoter.stop]
    if promoter.strand == '-':
        short_seq = str(long_seq)[promoter.stop:promoter.start]
        my_dna = Seq(short_seq)
        my_dna = my_dna.reverse_complement()
        short_seq=str(my_dna)
        if promoter.stop > len(my_dna)+promoter.start :
            lib.log.info("Promoter of gene " + promoter.locus_tag + " can't be fully extracted")
    
    short_seq_record = SeqRecord(Seq(short_seq), id=promoter.locus_tag, name=promoter.locus_tag, description=promoter.source+'_'+str(promoter.start)+':'+str(promoter.stop))
    return short_seq_record


"""
    Insert TFBS prediction
"""
def insert_tfbs_prediction(tfbs):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO tfbs VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                            (tfbs.regulatory_interaction_id, 
                                                            tfbs.pwm_id, 
                                                            tfbs.strand, 
                                                            tfbs.start, 
                                                            tfbs.end, 
                                                            tfbs.sequence, 
                                                            tfbs.weight, 
                                                            tfbs.pval, 
                                                            tfbs.ln_pval, 
                                                            tfbs.sig))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE tfbs")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE tfbs ")
        lib.log.info(str(tfbs.regulatory_interaction_id) + " " +
                    str(tfbs.pwm_id) + " " + 
                    str(tfbs.strand) + " " +
                    str(tfbs.sequence))
"""
    Construct GRN TFBS Predictions
"""
def construct_grn_tfbs_predictions(organism_id):
    tf_list = list()
    pwm_list = list()
    regulatory_interactions_list = list()
    
    # parse fasta file and turn into dictionary
    genome = SeqIO.to_dict(SeqIO.parse(open(in_file_genome), 'fasta'))
    
    # select TFs
    tf_list = select_tfs_by_organism(organism_id)
    for tf in tf_list:
        # Find PWMs
        pwm_list = select_pwms_by_locus_tag(tf.locus_tag)
        for pwm in pwm_list:
            #Find Regulatory Interactions
            regulatory_interactions_list = select_regulatory_interactions_by_tf_locus_tag(pwm.locus_tag)
            for regulatory_interaction in regulatory_interactions_list:
                
                prediction_file = out_folder + "tfbs_predictions/" + organism_id + "/" +regulatory_interaction.tf_locus_tag+"-"+regulatory_interaction.tg_locus_tag+"-"+pwm.motif_id+ ".txt"
                # this condition verify if exists a TFBS prediction already performed for this pwm
                # to avoid duplicated predictions with the same tg_promoter and pwm
                # note that pwms in Cis-Bp could be inferred by similarity tf_status = I
                if not os.path.isfile(prediction_file):
                    #Extract Promoter
                    promoter_sequence = extract_promoter(genome, regulatory_interaction.tg_locus_tag)
                    
                    ## Prepare matrix
                    in_matrix = in_folder + "Cis-BP/pwms/" + pwm.motif_id + ".txt"

                    # Input matrix (file content)
                    in_matrix_file = open(in_matrix, "r")
                    matrix = in_matrix_file.read()
                    in_matrix_file.close()

                    #Perform TFBS Prediction
                    ## Perform SOAP request on RSAT server with the matrix and the FASTA sequence
                    lib.log.info("Call RSAT Web Service: "+'\n'+
                             " tf_locus_tag: "+ regulatory_interaction.tf_locus_tag+'\n'+
                             " tg_locus_tag: "+ regulatory_interaction.tg_locus_tag+'\n'+
                             " pwm_id: "+ pwm.motif_id)

                    # Call web service matrix_scan
                    if matrix == ("Pos	A	C	G	T"+'\n'):
                        lib.log.info("Null Matrix: " + pwm.motif_id)
                    else:
                        result = call_matrix_scan(rsat_service, promoter_sequence, matrix)
                        # Write result in output file
                        with open(prediction_file, 'w') as out_file:
                            out_file.write(result)
                        out_file.close()

                        with open (prediction_file) as in_file:
                            for line in in_file:
                                if line.startswith("#"): continue
                                parts = line.strip().split("\t")
                                #Create new TFBS prediction for each RSAT prediction result
                                strand = urllib.parse.unquote(parts[3])
                                start = urllib.parse.unquote(parts[4])
                                end = urllib.parse.unquote(parts[5])
                                sequence = urllib.parse.unquote(parts[6])
                                weight = urllib.parse.unquote(parts[7])
                                pval = urllib.parse.unquote(parts[8])
                                ln_pval = urllib.parse.unquote(parts[9])
                                sig = urllib.parse.unquote(parts[10])
                                tfbs = TFBS(0, regulatory_interaction.id, pwm.id, strand, start, end, sequence, weight, pval, ln_pval, sig)
                                insert_tfbs_prediction(tfbs)
                            in_file.close()
                else:
                    lib.log.info("TFBS Prediction File already exists: " +regulatory_interaction.tf_locus_tag+"-"+regulatory_interaction.tg_locus_tag+"-"+pwm.motif_id+ ".txt")

"""
    Main function of this program
""" 
if __name__ == '__main__':
    
    #Create a database connection
    dbConnection = create_db_connection()
    
    #Load Input Files
    #parse_organism_file(in_file_organisms)
    #gff3_handler(in_file_genes)
    #parse_protein_file(in_file_proteins)
    #parse_pwm_file(in_file_pwms)
    #parse_model_regulatory_file(in_file_model_regulatory)
    
    #Construct GRN
    #####Orthology
    #execute ProteinOrtho *****
    #parse_orthology_file(in_file_orthology)
    #construct_grn_orthology('2', '4')
    #####TFBS Predictions
    construct_grn_tfbs_predictions('4')
    
    if (dbConnection):
        dbConnection.close()
        lib.log.info("The DB connection is closed")