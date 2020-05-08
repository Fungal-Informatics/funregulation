"""
###########################################################
#
# funregulation_network_inference.py
#
# FunRegulation: Network Inference (Alexandre Rafael Lenz)
# Universidade do Estado da Bahia (UNEB) / Universidade de Caxias do Sul (UCS)
# https://alexandrelenz@github.com/alexandrelenz/funregulation.git
# Last update: 07.05.2020 (Lenz, A. R.)
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
# # Python packages required:
# i)   pip3 install Biopython
# ii)  pip3 install suds.jurko
# iii) pip3 install natsort
# iv)  pip3 install pysqlite3
#
#
#   INPUT DATA:
#
# a) Genomic data: Four fungal genomes and proteomes were used:
#    - Penicillium ucsensis 2HH (this study);
#    - Penicillium oxalicum 114-2 (GCA000346795.1 pdev1.0 from GenBank);
#    - Aspergillus nidulansFGSC A4 (v. s10-m04-r06 from AspGD / GCF000149205.2 ASM14920v2 from GenBank);
#    - Neurospora crassaOR74A (v. 12 from Broad88Institute / GCF000182925.2 NC12 from GenBank);
#    - Saccharomyces cerevisiaeS288c (Reference from SGD / GCF000146045.2 R64 from GenBank).
#
# b) Regulatory interactions from related species:
#    Previous regulatory interactions of A. nidulans and N. crassa available in (Hu et al. 2018) 
#    and regulatory interactions of S. cerevisiae available in YEASTRACT (Monteiro et al. 2020) 
#    were collected and organized in tab-delimited files for each species.
#    # Regulatory interaction files:
#    i)   in_file_regulatory_Ani
#    ii)  in_file_regulatory_Ncr
#    iii) in_file_regulatory_Sce (6 partitioned files due the data volume)
#
# c) A cross-validation was performed to check locus tag and gene name for each regulatory interaction (b), 
#    crossing information from the reference genomes (GenBank) and regulatory interactions.
#    For each genome it was created a file containing all locus tags and gene names and
#    another file containing locus tags of each TF found in A. nidulans, N. crassa and S. cerevisiae.
#    # Gene names files:
#    i)   in_file_names_Puc
#    ii)  in_file_names_Pox
#    iii) in_file_names_Ani
#    iv)  in_file_names_Ncr
#    v)   in_file_names_Sce
#    # tf_locus_tag files:
#    vi)  in_file_tfs_Ani
#    vii) in_file_tfs_Ncr
#    viii)in_file_tfs_Sce
#
# d) ProteinOrtho (V6.0.15) Lechner et al. (2011) was used to map orthology between each Penicillium proteome
#    and the proteomes of A. nidulans, N. crassa and S. cerevisiae, individually.
#    # ProteinOrtho orthology files:
#    i)   in_file_ortho_Pox_Ani
#    ii)  in_file_ortho_Pox_Ncr
#    iii) in_file_ortho_Pox_Sce
#    iv)  in_file_ortho_Puc_Pox
#    v)   in_file_ortho_Puc_Ani
#    v)   in_file_ortho_Puc_Ncr
#    v)   in_file_ortho_Puc_Sce
#
# e) Information of Transcription Factors and PWMs were obtained from CIS-BP Database (Hu et al. 2013).
#    # CIS-BP information files:
#    i)   in_file_cisbp_Pox
#    ii)  in_file_cisbp_Ani
#    iii) in_file_cisbp_Ncr
#    iv)  in_file_cisbp_Sce
#   
# f) PWMs in 'cis-bp' format from P.oxalicum, A. nidulans, N. crassa and S. cerevisiae 
#    were obtained from CIS-BP Database (Hu et al. 2013).
#    # PWM files:
#    i) in_folder/pwms/
#   
# g) Regulatory sequences of P. ucsensis 2HH and P. oxalicum 114-2 were obtained by 
#    funregulation_promoter_extract.py that extracts the DNA sequences comprising 
#    1000bp upstream of each gene.
#    # Promoter regions in fasta format:
#    i)  in_folder/promoters/Poxalicum/
#    ii) in_folder/promoters/Pucsensis/
#   
#   SQLITE3 DATABASE:
#
#    The large volume and interconnectivity of data makes it difficult to access 
#    specific records in text files and makes it impossible to load complete files 
#    on machines with low memory capacity. Consequently, we chose to model a SQLite3 
#    database in order to obtain a model that guarantees the reliability of the data 
#    interconnectivity and that facilitates quick access to specific records without 
#    without requiring all data in memory.
#    # SQLite3 information files:
#    i) sqlite_db_folder
#    ii) sqlite_db_file_name
#    iii) sqlite_db_script_name
#          
#   OUTPUT FILES:
#
#    # Tab-delimited final output files:
#    i)   out_file_pucsensis_network
#    ii)  out_file_pucsensis_tfbs
#    iii) out_folder/tfbs_predictions/Pucsensis/
#    iv)  out_file_poxalicum_network
#    i)   out_file_poxalicum_tfbs
#    vi)  out_folder/tfbs_predictions/Poxalicum/
#
###########################################################
"""

import sys
import os, platform
import urllib.parse
import lib.library as lib
from collections import namedtuple
from suds.client import Client
from Bio import SeqIO
import sqlite3
import csv
from sqlite3 import Error

networkInferenceVersion = '0.3 - 02/05/2020'

"""
#################### INPUT PATHS ######################
"""
#Initialize Folder Paths
in_folder = "/Users/arlenz/NetBeansProjects/funregulation/GRN/input/"
out_folder = "/Users/arlenz/NetBeansProjects/funregulation/GRN/output/"

#Initialize Gene names files Paths
in_file_names_Puc = os.path.join(in_folder,'gene_names/Pucsensis_gene_names.tsv')
in_file_names_Pox = os.path.join(in_folder,'gene_names/Poxalicum_gene_names.tsv')
in_file_names_Ani = os.path.join(in_folder,'gene_names/Anidulans_gene_names.tsv')
in_file_names_Ncr = os.path.join(in_folder,'gene_names/Ncrassa_gene_names.tsv')
in_file_names_Sce = os.path.join(in_folder,'gene_names/Scerevisiae_gene_names.tsv')
#Initialize tf_locus_tag files Paths
in_file_tfs_Ani = os.path.join(in_folder,'tfs/Anidulans_tfs.tsv')
in_file_tfs_Ncr = os.path.join(in_folder,'tfs/Ncrassa_tfs.tsv')
in_file_tfs_Sce = os.path.join(in_folder,'tfs/Scerevisiae_tfs.tsv')
#Initialize Poxalicum proteinortho files Paths
in_file_ortho_Pox_Ani = os.path.join(in_folder,'orthology/Poxalicum/Poxalicum_Anidulans.proteinortho.tsv')
in_file_ortho_Pox_Ncr = os.path.join(in_folder,'orthology/Poxalicum/Poxalicum_Ncrassa.proteinortho.tsv')
in_file_ortho_Pox_Sce = os.path.join(in_folder,'orthology/Poxalicum/Poxalicum_Scerevisiae.proteinortho.tsv')
#Initialize Pucsensis proteinortho files Paths
in_file_ortho_Puc_Pox = os.path.join(in_folder,'orthology/Pucsensis/Pucsensis_Poxalicum.proteinortho.tsv')
in_file_ortho_Puc_Ani = os.path.join(in_folder,'orthology/Pucsensis/Pucsensis_Anidulans.proteinortho.tsv')
in_file_ortho_Puc_Ncr = os.path.join(in_folder,'orthology/Pucsensis/Pucsensis_Ncrassa.proteinortho.tsv')
in_file_ortho_Puc_Sce = os.path.join(in_folder,'orthology/Pucsensis/Pucsensis_Scerevisiae.proteinortho.tsv')
#Initialize CIS-BP information files Paths
in_file_cisbp_Pox = os.path.join(in_folder,'cis-bp/Poxalicum_cis-bp_information.tsv')
in_file_cisbp_Ani = os.path.join(in_folder,'cis-bp/Anidulans_cis-bp_information.tsv')
in_file_cisbp_Ncr = os.path.join(in_folder,'cis-bp/Ncrassa_cis-bp_information.tsv')
in_file_cisbp_Sce = os.path.join(in_folder,'cis-bp/Scerevisiae_cis-bp_information.tsv')
#Initialize regulatory interaction files Paths
in_file_regulatory_Ani = os.path.join(in_folder,'regulatory_interactions/Anidulans_regulatory_interactions.tsv')
in_file_regulatory_Ncr = os.path.join(in_folder,'regulatory_interactions/Ncrassa_regulatory_interactions.tsv')
in_file_regulatory_Sce_part1 = os.path.join(in_folder,'regulatory_interactions/Scerevisiae_regulatory_interactions_partition1.tsv')
in_file_regulatory_Sce_part2 = os.path.join(in_folder,'regulatory_interactions/Scerevisiae_regulatory_interactions_partition2.tsv')
in_file_regulatory_Sce_part3 = os.path.join(in_folder,'regulatory_interactions/Scerevisiae_regulatory_interactions_partition3.tsv')
in_file_regulatory_Sce_part4 = os.path.join(in_folder,'regulatory_interactions/Scerevisiae_regulatory_interactions_partition4.tsv')
in_file_regulatory_Sce_part5 = os.path.join(in_folder,'regulatory_interactions/Scerevisiae_regulatory_interactions_partition5.tsv')
in_file_regulatory_Sce_part6 = os.path.join(in_folder,'regulatory_interactions/Scerevisiae_regulatory_interactions_partition6.tsv')

#Initialize Output File Paths
out_file_pucsensis_network = os.path.join(out_folder,'network/Pucsensis_network.tsv')
out_file_pucsensis_tfbs = os.path.join(out_folder,'tfbs_predictions/Pucsensis_tfbs_prediction.tsv')
out_file_poxalicum_network = os.path.join(out_folder,'network/Poxalicum_network.tsv')
out_file_poxalicum_tfbs = os.path.join(out_folder,'tfbs_predictions/Poxalicum_tfbs_prediction.tsv')

#SQLite3 definitions
sqlite_db_folder = '/Users/arlenz/NetBeansProjects/funregulation/GRN/'
sqlite_db_file_name = 'GRN.sqlite3'
sqlite_db_file = os.path.join(sqlite_db_folder,sqlite_db_file_name)
sqlite_db_script_name = 'funregulation_GRN.sql'
sqlite_db_script_file = os.path.join(sqlite_db_folder,sqlite_db_script_name)
sqliteConnection = None;

"""
#################### PRE-PROCESSING ######################
"""
#create folder structure
if not os.path.isdir(out_folder):
    os.makedirs(out_folder)
    os.makedirs(os.path.join(out_folder, 'tfbs_predictions'))
    os.makedirs(os.path.join(out_folder, 'logfiles'))
else:
    #make sure subdirectories exist
    dirs = [os.path.join(out_folder, 'tfbs_predictions'), os.path.join(out_folder, 'logfiles')]
    for d in dirs:
        if not os.path.isdir(d):
            os.makedirs(d)

#create log file
log_name = os.path.join(out_folder, 'logfiles', 'funregulation-networkinference.log')
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

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
wsdlUrl = 'http://embnet.ccg.unam.mx/rsat/web_services/RSATWS2.wsdl'

# Create the client
try:
    client = Client(wsdlUrl)
except Error as e:
    lib.log.info("Connection to RSAT Server failed!", e)
    raise
# Need the service interface to perform requests
rsat_service = client.service

"""
    Define client header (optional)
"""
userAgent = 'RSAT-Client/v%s (%s; Python %s; %s)' % (
    networkInferenceVersion, 
    os.path.basename( __file__ ),
    platform.python_version(), 
    platform.system()
)
httpHeaders = {'User-agent': userAgent}
client.set_options(headers=httpHeaders)
client.set_options(timeout=300)

#RSAT config parameters
uth_pval = "0.0001"
matrix_format = "cis-bp"

"""
#######################   CLASSES    ####################
"""
"""
    Class defining a Gene
"""
class Gene:
    def __init__(self, organism, locus_tag, name, is_tf):
        self.organism = organism
        self.locus_tag = locus_tag
        self.name = name
        self.is_tf = is_tf

"""
    Class defining a Ortholog relationship
"""
class Ortho:
    def __init__(self, src_gene, ortho_gene):
        self.src_gene = src_gene
        self.ortho_gene = ortho_gene

"""
    Class defining a PWM matrix from CIS-BP
"""
class Pwm:
    def __init__(self, motif_id, tf_status, family_name, motif_type, msource_author, msource_year, pmid):
        self.motif_id = motif_id
        self.tf_status = tf_status
        self.family_name = family_name
        self.motif_type = motif_type
        self.msource_author = msource_author
        self.msource_year = msource_year
        self.pmid = pmid

"""
    Class defining a Regulatory Interaction
"""
class RegulatoryInteraction:
    def __init__(self, tf_ortho, tg_ortho, interaction):
        self.tf_ortho = tf_ortho
        self.tg_ortho = tg_ortho
        self.interaction = interaction
        self.tfbs_list = list()
    
    def add_tfbs_prediction(self, tfbs_prediction):
        self.tfbs_list.append(tfbs_prediction)

"""
    Class defining a TFBS Prediction
"""
class TFBSPrediction:
    def __init__(self, pwm, strand, start, end, sequence, weight, pval, ln_pval, sig):
        self.pwm = pwm
        self.strand = strand
        self.start = start
        self.end = end
        self.sequence = sequence
        self.weight = weight
        self.pval = pval
        self.ln_pval = ln_pval
        self.sig = sig

"""
    Class defining a Global Gene Regulatory Network Node
"""
class NetworkNode:
    def __init__(self, tf_gene, tf_ortho_names, tg_gene, tg_ortho_names, interaction, ortho_pox, ortho_ani, ortho_ncr, ortho_sce, tfbs_count):
        self.tf_gene = tf_gene
        self.tf_ortho_names = tf_ortho_names
        self.tg_gene = tg_gene
        self.tg_ortho_names = tg_ortho_names
        self.interaction = interaction
        self.ortho_pox = ortho_pox
        self.ortho_ani = ortho_ani
        self.ortho_ncr = ortho_ncr
        self.ortho_sce = ortho_sce
        self.tfbs_count = tfbs_count

"""
#################### SQLite3 DEFINITIONS ###################### 
"""

"""
    SQLite3 Delete Current Database and Create a new empty Database
"""
def create_sqlite_database(sqlite_db_file,sqlite_db_script_file):
    # Delete the old database
    if os.path.isfile(sqlite_db_file):
        try:
            os.remove(sqlite_db_file)
            lib.log.info("SQLite3 Database Successfully deleted: " + sqlite_db_file)
        except Exception as e:
            lib.log.info(errorMessage = sqlite_db_file + ': ' + str(e))
            raise
    #Create new the database
    qry = open(sqlite_db_script_file, 'r').read()
    sqlite3.complete_statement(qry)
    conn = sqlite3.connect(sqlite_db_file)
    cursor = conn.cursor()
    try:
        cursor.executescript(qry)
        lib.log.info("SQLite3 Database Successfully created: " + sqlite_db_file)
    except Exception as e:
        lib.log.info(errorMessage = sqlite_db_file + ': ' + str(e))
        cursor.close()
        raise

""" create a database connection to a SQLite database
"""
def create_sqlite_connection(sqlite_db_file):
    try:
        conn = sqlite3.connect(sqlite_db_file)
        lib.log.info("Successfully Connected to SQLite: " + sqlite3.version)
        return conn
    except Error as e:
        lib.log.info(e)

"""
#################### SQLite3 INSERTS ###################### 
"""

"""
    SQLite3 Insert gene
"""
def sqlite_insert_gene(gene):
    try:
        cursor = sqliteConnection.cursor()
        sqlite_insert_query = "INSERT INTO gene (organism,locus_tag,name,is_tf) VALUES(?, ?, ?, ?)"
        count = cursor.execute(sqlite_insert_query,(str(gene.organism), 
                                                    str(gene.locus_tag),
                                                    str(gene.name),
                                                    int(gene.is_tf)))
        sqliteConnection.commit()
        lib.log.info("Record inserted successfully into "+ sqlite_db_file_name + " TABLE gene ")
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to insert data into sqlite table", error)
        lib.log.info(str(gene.organism) + " " +
                    str(gene.locus_tag) + " " + 
                    str(gene.name) + " " +
                    str(gene.is_tf))

"""
    SQLite3 Insert ortho
"""
def sqlite_insert_ortho(ortho):
    try:
        cursor = sqliteConnection.cursor()
        sqlite_insert_query = "INSERT INTO ortho (src_organism,src_locus_tag,ortho_organism,ortho_locus_tag) VALUES(?, ?, ?, ?)"
        count = cursor.execute(sqlite_insert_query,(str(ortho.src_gene.organism), 
                                                    str(ortho.src_gene.locus_tag),
                                                    str(ortho.ortho_gene.organism),
                                                    str(ortho.ortho_gene.locus_tag)))
        sqliteConnection.commit()
        lib.log.info("Record inserted successfully into "+ sqlite_db_file_name + " TABLE ortho ")
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to insert data into sqlite table", error)
        lib.log.info(str(ortho.src_gene.organism) + " " +
                    str(ortho.src_gene.locus_tag) + " " + 
                    str(ortho.ortho_gene.organism) + " " +
                    str(ortho.ortho_gene.locus_tag))
        
"""
    SQLite3 Insert pwm
"""
def sqlite_insert_pwm(gene, pwm):
    try:
        cursor = sqliteConnection.cursor()
        sqlite_insert_query = "INSERT INTO pwm (organism,locus_tag,motif_id,tf_status,family_name,motif_type,msource_author,msource_year,pmid) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)"
        count = cursor.execute(sqlite_insert_query,(str(gene.organism), 
                                                    str(gene.locus_tag),
                                                    str(pwm.motif_id),
                                                    str(pwm.tf_status),
                                                    str(pwm.family_name),
                                                    str(pwm.motif_type),
                                                    str(pwm.msource_author),
                                                    str(pwm.msource_year),
                                                    str(pwm.pmid)))
        sqliteConnection.commit()
        lib.log.info("Record inserted successfully into "+ sqlite_db_file_name + " TABLE pwm ")
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to insert data into sqlite table", error)
        lib.log.info(str(gene.organism) + " " +
                    str(gene.locus_tag) + " " +
                    str(pwm.motif_id) + " " +
                    str(pwm.tf_status) + " " +
                    str(pwm.family_name) + " " +
                    str(pwm.motif_type) + " " +
                    str(pwm.msource_author) + " " +
                    str(pwm.msource_year) + " " +
                    str(pwm.pmid))
        
"""
    SQLite3 Insert regulation
"""
def sqlite_insert_regulatory_interaction(interaction):
    try:
        cursor = sqliteConnection.cursor()
        sqlite_insert_query = "INSERT INTO regulation (src_organism,src_tf_locus_tag,src_tg_locus_tag,interaction,ortho_organism,ortho_tf_locus_tag,ortho_tg_locus_tag) VALUES(?, ?, ?, ?, ?, ?, ?)"
        count = cursor.execute(sqlite_insert_query,(str(interaction.tf_ortho.src_gene.organism), 
                                                    str(interaction.tf_ortho.src_gene.locus_tag),
                                                    str(interaction.tg_ortho.src_gene.locus_tag),
                                                    str(interaction.interaction),
                                                    str(interaction.tf_ortho.ortho_gene.organism),
                                                    str(interaction.tf_ortho.ortho_gene.locus_tag),
                                                    str(interaction.tg_ortho.ortho_gene.locus_tag)))
        sqliteConnection.commit()
        lib.log.info("Record inserted successfully into "+ sqlite_db_file_name + " TABLE regulation ")
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to insert data into sqlite table", error)
        lib.log.info(str(interaction.tf_ortho.src_gene.organism) + " " +
                    str(interaction.tf_ortho.src_gene.locus_tag) + " " +
                    str(interaction.tg_ortho.src_gene.locus_tag) + " " +
                    str(interaction.interaction) + " " +
                    str(interaction.tf_ortho.ortho_gene.organism) + " " +
                    str(interaction.tf_ortho.ortho_gene.locus_tag) + " " +
                    str(interaction.tg_ortho.ortho_gene.locus_tag))
                    
"""
    SQLite3 Insert regulation
"""
def sqlite_insert_tfbs_prediction(interaction, tfbs):
    try:
        cursor = sqliteConnection.cursor()
        sqlite_insert_query = "INSERT INTO tfbs_prediction (src_organism,src_tf_locus_tag,src_tg_locus_tag,ortho_organism,ortho_tf_locus_tag,motif_id,strand,start,end,sequence,weight,pval,ln_pval,sig) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
        count = cursor.execute(sqlite_insert_query,(str(interaction.tf_ortho.src_gene.organism), 
                                                    str(interaction.tf_ortho.src_gene.locus_tag),
                                                    str(interaction.tg_ortho.src_gene.locus_tag),
                                                    str(interaction.tf_ortho.ortho_gene.organism),
                                                    str(interaction.tf_ortho.ortho_gene.locus_tag),
                                                    str(tfbs.pwm),
                                                    str(tfbs.strand),
                                                    int(tfbs.start),
                                                    int(tfbs.end),
                                                    str(tfbs.sequence),
                                                    float(tfbs.weight),
                                                    str(tfbs.pval),
                                                    float(tfbs.ln_pval),
                                                    float(tfbs.sig)))
        sqliteConnection.commit()
        lib.log.info("Record inserted successfully into "+ sqlite_db_file_name + " TABLE tfbs_prediction ")
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to insert data into sqlite table", error)
        lib.log.info(str(interaction.tf_ortho.src_gene.organism) + " " +
                    str(interaction.tf_ortho.src_gene.locus_tag) + " " +
                    str(interaction.tf_ortho.ortho_gene.organism) + " " +
                    str(interaction.tf_ortho.ortho_gene.locus_tag) + " " +
                    str(interaction.tg_ortho.ortho_gene.locus_tag) + " " +
                    str(tfbs.pwm.motif_id) + " " +
                    str(tfbs.strand) + " " +
                    str(tfbs.start) + " " +
                    str(tfbs.end) + " " +
                    str(tfbs.sequence) + " " +
                    str(tfbs.weight) + " " +
                    str(tfbs.pval) + " " +
                    str(tfbs.ln_pval) + " " +
                    str(tfbs.sig))

"""
    SQLite3 Insert Network Node
"""
def sqlite_insert_network_node(interaction):
    tf_ortho_names = ""
    tg_ortho_names = ""
    if interaction.tf_ortho.ortho_gene.name != interaction.tf_ortho.ortho_gene.locus_tag:
        tf_ortho_names = interaction.tf_ortho.ortho_gene.name
    if interaction.tg_ortho.ortho_gene.name != interaction.tg_ortho.ortho_gene.locus_tag:
        tg_ortho_names = interaction.tg_ortho.ortho_gene.name

    ortho_pox = 0
    ortho_ani = 0
    ortho_ncr = 0
    ortho_sce = 0
    if interaction.tf_ortho.ortho_gene.organism == "Poxalicum":
        ortho_pox = 1
    elif interaction.tf_ortho.ortho_gene.organism == "Anidulans":
        ortho_ani = 1
    elif interaction.tf_ortho.ortho_gene.organism == "Ncrassa":
        ortho_ncr = 1
    elif interaction.tf_ortho.ortho_gene.organism == "Scerevisiae":
        ortho_sce = 1
        
    tfbs_count = get_tfbs_count(interaction.tf_ortho.src_gene.locus_tag,interaction.tg_ortho.src_gene.locus_tag)
    
    try:
        cursor = sqliteConnection.cursor()
        sqlite_insert_query = "INSERT INTO network_node (organism,tf_locus_tag,tf_name,tf_ortho_names,tg_locus_tag,tg_name,tg_ortho_names,interaction,ortho_pox,ortho_ani,ortho_ncr,ortho_sce,tfbs_count) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
        count = cursor.execute(sqlite_insert_query,(str(interaction.tf_ortho.src_gene.organism),
                                                    str(interaction.tf_ortho.src_gene.locus_tag),
                                                    str(interaction.tf_ortho.src_gene.name),
                                                    str(tf_ortho_names),
                                                    str(interaction.tg_ortho.src_gene.locus_tag),
                                                    str(interaction.tg_ortho.src_gene.name),
                                                    str(tg_ortho_names),
                                                    str(interaction.interaction),
                                                    int(ortho_pox),
                                                    int(ortho_ani),
                                                    int(ortho_ncr),
                                                    int(ortho_sce),
                                                    int(tfbs_count)))
        sqliteConnection.commit()
        lib.log.info("Record inserted successfully into "+ sqlite_db_file_name + " TABLE network_node ")
        lib.log.info(str(interaction.tf_ortho.src_gene.locus_tag) +"  "+ str(interaction.tg_ortho.src_gene.locus_tag))
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to insert data into sqlite table", error)
        lib.log.info(str(node.tf_gene.organism) + " " +
                    str(node.tf_gene.locus_tag) + " " +
                    str(node.tf_gene.name) + " " +
                    str(node.tg_gene.locus_tag) + " " +
                    str(node.tg_gene.name) + " " +
                    str(node.interaction) + " " +
                    str(node.ortho_pox) + " " +
                    str(node.ortho_ani) + " " +
                    str(node.ortho_ncr) + " " +
                    str(node.ortho_sce) + " " +
                    str(node.tfbs_count))

"""
    ####################### SQLite3 UPDATES #######################
"""

"""
    SQLite3 Update TF Gene 
"""
def sqlite_update_gene_tf(locus_tag,is_tf):
    try:
        cursor = sqliteConnection.cursor()
        sqlite_update_query = "UPDATE gene SET is_tf = ? WHERE locus_tag = ?"
        count = cursor.execute(sqlite_update_query,(is_tf, locus_tag))
        sqliteConnection.commit()
        lib.log.info("Record updated successfully into "+ sqlite_db_file_name + " TABLE gene ")
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to update data into sqlite table", error)
        lib.log.info(locus_tag + " " + str(gene.is_tf))

"""
    SQLite3 Update Network Node 
"""
def sqlite_update_network_node(node,interaction):
    ortho_pox = node.ortho_pox
    ortho_ani = node.ortho_ani
    ortho_ncr = node.ortho_ncr
    ortho_sce = node.ortho_sce
    if interaction.tf_ortho.ortho_gene.organism == "Poxalicum":
        ortho_pox = 1
    elif interaction.tf_ortho.ortho_gene.organism == "Anidulans":
        ortho_ani = 1
    elif interaction.tf_ortho.ortho_gene.organism == "Ncrassa":
        ortho_ncr = 1
    elif interaction.tf_ortho.ortho_gene.organism == "Scerevisiae":
        ortho_sce = 1
    
    tf_ortho_names = node.tf_ortho_names
    if interaction.tf_ortho.ortho_gene.name != interaction.tf_ortho.ortho_gene.locus_tag and (interaction.tf_ortho.ortho_gene.name not in tf_ortho_names):
        tf_ortho_names = tf_ortho_names +','+ interaction.tf_ortho.ortho_gene.name
        
    tg_ortho_names = node.tg_ortho_names
    if (interaction.tg_ortho.ortho_gene.name != interaction.tg_ortho.ortho_gene.locus_tag) and (interaction.tg_ortho.ortho_gene.name not in tg_ortho_names):
        tg_ortho_names = tg_ortho_names +','+ interaction.tg_ortho.ortho_gene.name
    
    tfbs_count = get_tfbs_count(node.tf_gene.locus_tag,node.tg_gene.locus_tag)
    
    try:
        cursor = sqliteConnection.cursor()
        sqlite_update_query = "UPDATE network_node SET tf_ortho_names = ? , tg_ortho_names = ? , ortho_pox = ? , ortho_ani = ? , ortho_ncr = ? , ortho_sce = ? , tfbs_count = ? WHERE tf_locus_tag = ? AND tg_locus_tag = ?"
        count = cursor.execute(sqlite_update_query,(str(tf_ortho_names),
                                                    str(tg_ortho_names),
                                                    int(ortho_pox),
                                                    int(ortho_ani),
                                                    int(ortho_ncr),
                                                    int(ortho_sce),
                                                    int(tfbs_count),
                                                    str(node.tf_gene.locus_tag),
                                                    str(node.tg_gene.locus_tag)))
        sqliteConnection.commit()
        lib.log.info("Record updated successfully into "+ sqlite_db_file_name + " TABLE network_node ")
        lib.log.info(str(node.tf_gene.locus_tag) +"  "+ str(node.tg_gene.locus_tag))
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to update data into sqlite table", error)
        lib.log.info(str(node.tf_gene.locus_tag)+ " " + 
                    str(node.tg_gene.locus_tag) + " " +
                    str(tf_ortho_names)+ " " + 
                    str(tg_ortho_names)+ " " + 
                    str(ortho_pox)+ " " + 
                    str(ortho_ani)+ " " + 
                    str(ortho_ncr)+ " " + 
                    str(ortho_sce)+ " " + 
                    str(tfbs_count))

"""
    Update network nodes of Pucsensis x ortho Poxalicum
"""
def sqlite_update_network_node_puc_by_pox_ortho():
    data_nodes_puc = None
    data_tf_ortho = None
    data_tg_ortho = None
    try:
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM network_node WHERE organism = ?",[str("Pucsensis")])
        data_nodes_puc = cursor.fetchall()
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to select data into sqlite table", error)
        
    for row in data_nodes_puc:
        try:
            cursor = sqliteConnection.cursor()
            cursor.execute("SELECT * FROM ortho WHERE src_organism = ? AND src_locus_tag = ? AND ortho_organism = ?",[str("Pucsensis"),str(row[1]),str("Poxalicum")])
            data_tf_ortho = cursor.fetchall()
            cursor.execute("SELECT * FROM ortho WHERE src_organism = ? AND src_locus_tag = ? AND ortho_organism = ?",[str("Pucsensis"),str(row[4]),str("Poxalicum")])
            data_tg_ortho = cursor.fetchall()
            cursor.close()
        except sqlite3.Error as error:
            lib.log.info("Failed to select data into sqlite table", error)
            
        if len(data_tf_ortho) !=0 and len(data_tg_ortho) !=0 :
            try:
                cursor = sqliteConnection.cursor()
                sqlite_update_query = "UPDATE network_node SET ortho_pox = ? WHERE organism = ? AND tf_locus_tag = ? AND tg_locus_tag = ?"
                count = cursor.execute(sqlite_update_query,(int(1),
                                                            str("Pucsensis"),
                                                            str(str(row[1])),
                                                            str(str(row[4]))))
                sqliteConnection.commit()
                lib.log.info("Record updated successfully into "+ sqlite_db_file_name + " TABLE network_node ")
                lib.log.info(str(row[1]) +"  "+ str(row[4]))
                cursor.close()
            except sqlite3.Error as error:
                lib.log.info("Failed to insert data into sqlite table", error)
                lib.log.info("Pucsensis" + " " +
                            str(row[1]) + " " +
                            str(row[4]))

"""
    ####################### SQLite3 SELECTS #######################
"""
"""
    SQLite3 Select Orthologs by src_organism e ortho_locus_tag return list
"""
def sqlite_get_orthologs_by_src_organism(src_organism, ortho_locus_tag):
    try:
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM ortho WHERE src_organism = ? AND ortho_locus_tag = ?",[str(src_organism), str(ortho_locus_tag)])
        records=cursor.fetchall()
        ortho_list = list()
        if len(records)==0:
            return ortho_list
        else:
            for row in records:
                
                src_gene = sqlite_get_gene_by_locus_tag(row[1])
                ortho_gene = sqlite_get_gene_by_locus_tag(row[3])
                
                ortholog = Ortho(src_gene,ortho_gene)
                
                ortho_list.append(ortholog)
            return ortho_list
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to execute the select into sqlite table", error)
        lib.log.info(locus_tag)

"""
    SQLite3 Select Orthologs by ortho_organism e src_locus_tag return list
"""
def sqlite_get_orthologs_by_ortho_organism(ortho_organism, src_locus_tag):
    try:
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM ortho WHERE ortho_organism = ? AND src_locus_tag = ?",[str(ortho_organism), str(src_locus_tag)])
        records=cursor.fetchall()
        ortho_list = list()
        if len(records)==0:
            return ortho_list
        else:
            for row in records:
                
                src_gene = sqlite_get_gene_by_locus_tag(row[1])
                ortho_gene = sqlite_get_gene_by_locus_tag(row[3])
                
                ortholog = Ortho(src_gene,ortho_gene)
                
                ortho_list.append(ortholog)
            return ortho_list
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to execute the select into sqlite table", error)
        lib.log.info(locus_tag)
        
"""
    SQLite3 Select Gene by locus_tag
"""
def sqlite_get_gene_by_locus_tag(locus_tag):
    try:
        cursor = sqliteConnection.cursor()
        organism = cursor.execute("SELECT organism FROM gene WHERE locus_tag = ?", [str(locus_tag)]).fetchone()[0]
        name = cursor.execute("SELECT name FROM gene WHERE locus_tag = ?", [str(locus_tag)]).fetchone()[0]
        is_tf = cursor.execute("SELECT is_tf FROM gene WHERE locus_tag = ?", [str(locus_tag)]).fetchone()[0]
        gene = Gene(organism,locus_tag,name,is_tf)
        return gene
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to execute the select into sqlite table", error)
        lib.log.info(locus_tag)

"""
    SQLite3 Select genes by organism and is_tf
"""
def sqlite_get_tfs_by_organism(organism):
    try:
        cursor = sqliteConnection.cursor()
        records = cursor.execute("SELECT * FROM gene WHERE organism = ? AND is_tf = 1", [str(organism)]).fetchall()
        gene_list = list()
        for row in records:
            gene = Gene(row[0],row[1],row[2],row[3])
            gene_list.append(gene)
        lib.log.info("TFs by organism selected successfully")
        return gene_list
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to execute the select into sqlite table", error)
        lib.log.info(organism)
        
"""
    SQLite3 Select pwm by tf_locus_tag
"""
def sqlite_get_pwm_by_tf_locus_tag(tf_locus_tag):
    try:
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM pwm WHERE locus_tag = ?",[str(tf_locus_tag)])
        records=cursor.fetchall()
        pwm_list = list()
        for row in records:
            pwm = Pwm(row[2],row[3],row[4],row[5],row[6],row[7],row[8])
            pwm_list.append(pwm)
        return pwm_list
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to execute the select into sqlite table", error)
        lib.log.info(locus_tag + " " + motif_id)

"""
    SQLite3 Select Regulatory Interactions by tf_locus_tag
"""
def sqlite_get_interactions_by_tf_locus_tag(tf_locus_tag):
    try:
        lib.log.info("Looking for Interactions by tf_locus_tag: "+tf_locus_tag)
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM regulation WHERE src_tf_locus_tag = ?",[str(tf_locus_tag)])
        data = cursor.fetchall()
        interactions_list = list()
        for row in data:
            src_tf = sqlite_get_gene_by_locus_tag(row[1])
            ortho_tf = sqlite_get_gene_by_locus_tag(row[5])
            ortholog_tf = Ortho(src_tf,ortho_tf)
            
            src_tg = sqlite_get_gene_by_locus_tag(row[2])
            ortho_tg = sqlite_get_gene_by_locus_tag(row[6])
            ortholog_tg = Ortho(src_tg,ortho_tg)
            
            interaction = RegulatoryInteraction(ortholog_tf, ortholog_tg, row[3])
            interactions_list.append(interaction)
        lib.log.info("Interactions by tf_locus_tag selected successfully")
        return interactions_list
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to select data into sqlite table", error)

"""
    SQLite3 Select Regulatory Interactions by src_organism
"""
def sqlite_get_interactions_by_src_organism(src_organism):
    try:
        lib.log.info("Looking for Interactions by src_organism. This can take a long time.")
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM regulation WHERE src_organism = ?",[str(src_organism)])
        data = cursor.fetchall()
        interactions_list = list()
        for row in data:
            src_tf = sqlite_get_gene_by_locus_tag(row[1])
            ortho_tf = sqlite_get_gene_by_locus_tag(row[5])
            ortholog_tf = Ortho(src_tf,ortho_tf)
            
            src_tg = sqlite_get_gene_by_locus_tag(row[2])
            ortho_tg = sqlite_get_gene_by_locus_tag(row[6])
            ortholog_tg = Ortho(src_tg,ortho_tg)
            
            interaction = RegulatoryInteraction(ortholog_tf, ortholog_tg, row[3])
            interactions_list.append(interaction)
        lib.log.info("Interactions by src_organism selected successfully")
        return interactions_list
        cursor.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to select data into sqlite table", error)

"""
    SQLite3 Get Network Node by tf_locus_tag and tg_locus_tag
"""
def sqlite_get_network_node_by_tf_tg(tf_locus_tag, tg_locus_tag):
    try:
        lib.log.info("Looking for Network Node by tf_locus_tag and tg_locus_tag.")
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM network_node WHERE tf_locus_tag = ? AND tg_locus_tag = ?",[str(tf_locus_tag),str(tg_locus_tag)])
        data = cursor.fetchall()
        
        if(len(data) == 0):
            return None
        else:
            for row in data:
                tf_gene = sqlite_get_gene_by_locus_tag(row[1])
                tf_ortho_names = row[3]
                tg_gene = sqlite_get_gene_by_locus_tag(row[4])
                tg_ortho_names = row[6]
                interaction = row[7]
                ortho_pox = row[8] 
                ortho_ani = row[9]
                ortho_ncr = row[10] 
                ortho_sce = row[11]
                tfbs_count = row[12]
                
                node = NetworkNode(tf_gene,tf_ortho_names,tg_gene,tg_ortho_names,interaction,ortho_pox,ortho_ani,ortho_ncr,ortho_sce,tfbs_count)
                lib.log.info("Network Node by tf_locus_tag and tg_locus_tag selected successfully")
                cursor.close()
                return node
    except sqlite3.Error as error:
        lib.log.info("Failed to select data into sqlite table", error)


def get_tfbs_count(tf_locus_tag,tg_locus_tag):
    try:
        lib.log.info("Counting tfbs predictions by tf_locus_tag and tg_locus_tag.")
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM tfbs_prediction WHERE src_tf_locus_tag = ? AND src_tg_locus_tag = ?",[str(tf_locus_tag),str(tg_locus_tag)])
        data = cursor.fetchall()
        lib.log.info(str(len(data)) + "tfbs predicted for tf_locus_tag and tg_locus_tag")
        cursor.close()
        return len(data)
    except sqlite3.Error as error:
        lib.log.info("Failed to select data into sqlite table", error)

"""
    ####################### INPUT FILES #######################
"""

"""
    Load gene names input file
"""
def parse_name_file(organism, filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            locus_tag = urllib.parse.unquote(parts[0])
            gene_name = urllib.parse.unquote(parts[1])
            gene = Gene(organism,locus_tag,gene_name,0)
            sqlite_insert_gene(gene)
            
    in_file.close()
    lib.log.info(filename + " parsed correctly")

"""
    Load tfs input file
"""
def parse_tfs_file(filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            locus_tag = urllib.parse.unquote(parts[0])
            sqlite_update_gene_tf(locus_tag,1)
    in_file.close()
    lib.log.info(filename + " parsed correctly")
    
"""
    Load orthology input file
"""
def parse_orthology_file(filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            line_parts = line.strip().split("\t")
            
            source = urllib.parse.unquote(line_parts[3])
            ortholog = urllib.parse.unquote(line_parts[4])
            
            src_parts = source.strip().split(",")
            ortho_parts = ortholog.strip().split(",")
            
            for record_src in src_parts:
                for record_ortho in ortho_parts:
                    if (record_src != '*' and record_ortho != '*'):
                        src_gene = sqlite_get_gene_by_locus_tag(record_src)
                        ortho_gene = sqlite_get_gene_by_locus_tag(record_ortho)
                        ortho = Ortho(src_gene,ortho_gene)
                        sqlite_insert_ortho(ortho)
                    
    in_file.close()
    lib.log.info(filename + " parsed correctly")
    
"""
    Load CIS-BP PWM information input file
"""
def parse_pwm_information_file(filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            line_parts = line.strip().split("\t")
            
            motif_id = urllib.parse.unquote(line_parts[3]) 
            if motif_id != '.':
                locus_tag = urllib.parse.unquote(line_parts[5])
                gene = sqlite_get_gene_by_locus_tag(locus_tag)
                tf_status = urllib.parse.unquote(line_parts[8])
                family_name = urllib.parse.unquote(line_parts[9])
                motif_type = urllib.parse.unquote(line_parts[14])
                msource_author = urllib.parse.unquote(line_parts[17])
                msource_year = urllib.parse.unquote(line_parts[18])
                pmid = urllib.parse.unquote(line_parts[19])
                pwm = Pwm(motif_id, tf_status, family_name, motif_type, msource_author, msource_year, pmid)
                sqlite_insert_pwm(gene,pwm)
    in_file.close()
    lib.log.info(filename + " parsed correctly")
    
"""
    Load Regulatory Interactions files
"""
def parse_regulatory_interactions_file(src_organism, filename):
    lib.log.info("Parsing " + filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            
            tf_locus_tag = urllib.parse.unquote(parts[0])
            
            #update gene as TF if it was not at the input list
            gene = sqlite_get_gene_by_locus_tag(tf_locus_tag)
            if gene.is_tf == 0 :
                sqlite_update_gene_tf(tf_locus_tag,1)
            
            tg_locus_tag = urllib.parse.unquote(parts[1])
            regulatory_function = urllib.parse.unquote(parts[2])
            
            tf_orthologs = sqlite_get_orthologs_by_src_organism(src_organism, tf_locus_tag)
            tg_orthologs = sqlite_get_orthologs_by_src_organism(src_organism, tg_locus_tag)

            if len(tf_orthologs)!=0 and len(tg_orthologs)!=0:
                for record_tf in tf_orthologs:
                    for record_tg in tg_orthologs:
                        regulation = RegulatoryInteraction(record_tf, record_tg, regulatory_function)
                        if record_tf.src_gene.is_tf == 0 :
                            sqlite_update_gene_tf(str(record_tf.src_gene.locus_tag),1)
                        sqlite_insert_regulatory_interaction(regulation)
    in_file.close()
    lib.log.info(filename + " parsed correctly")

"""
    Load all input files
"""

def load_input_files():
    #Load gene_names files
    parse_name_file('Pucsensis', in_file_names_Puc)
    parse_name_file('Poxalicum', in_file_names_Pox)
    parse_name_file('Anidulans', in_file_names_Ani)
    parse_name_file('Ncrassa', in_file_names_Ncr)
    parse_name_file('Scerevisiae', in_file_names_Sce)
    
    #Load tf_locus_tag files
    parse_tfs_file(in_file_tfs_Ani)
    parse_tfs_file(in_file_tfs_Ncr)
    parse_tfs_file(in_file_tfs_Sce)

    #Load orthology files
    parse_orthology_file(in_file_ortho_Pox_Ani)
    parse_orthology_file(in_file_ortho_Pox_Ncr)
    parse_orthology_file(in_file_ortho_Pox_Sce)
    parse_orthology_file(in_file_ortho_Puc_Pox)
    parse_orthology_file(in_file_ortho_Puc_Ani)
    parse_orthology_file(in_file_ortho_Puc_Ncr)
    parse_orthology_file(in_file_ortho_Puc_Sce)

    #Load CIS-BP PWM information files
    parse_pwm_information_file(in_file_cisbp_Pox)
    parse_pwm_information_file(in_file_cisbp_Ani)
    parse_pwm_information_file(in_file_cisbp_Ncr)
    parse_pwm_information_file(in_file_cisbp_Sce)
    
    #Load regulatory intaractions files
    parse_regulatory_interactions_file('Poxalicum', in_file_regulatory_Ani)
    parse_regulatory_interactions_file('Poxalicum', in_file_regulatory_Ncr)
    parse_regulatory_interactions_file('Poxalicum', in_file_regulatory_Sce_part1)
    parse_regulatory_interactions_file('Poxalicum', in_file_regulatory_Sce_part2)
    parse_regulatory_interactions_file('Poxalicum', in_file_regulatory_Sce_part3)
    parse_regulatory_interactions_file('Poxalicum', in_file_regulatory_Sce_part4)
    parse_regulatory_interactions_file('Poxalicum', in_file_regulatory_Sce_part5)
    parse_regulatory_interactions_file('Poxalicum', in_file_regulatory_Sce_part6)
    
    parse_regulatory_interactions_file('Pucsensis', in_file_regulatory_Ani)
    parse_regulatory_interactions_file('Pucsensis', in_file_regulatory_Ncr)
    parse_regulatory_interactions_file('Pucsensis', in_file_regulatory_Sce_part1)
    parse_regulatory_interactions_file('Pucsensis', in_file_regulatory_Sce_part2)
    parse_regulatory_interactions_file('Pucsensis', in_file_regulatory_Sce_part3)
    parse_regulatory_interactions_file('Pucsensis', in_file_regulatory_Sce_part4)
    parse_regulatory_interactions_file('Pucsensis', in_file_regulatory_Sce_part5)
    parse_regulatory_interactions_file('Pucsensis', in_file_regulatory_Sce_part6)
    
"""
    ####################### PROCESSING #######################
"""
"""
    Update is_tf in Poxalicum and Pucsensis by orthology from known TFs
"""
def update_tfs():
    # Anidulans tf orthologs
    tf_genes_ani = sqlite_get_tfs_by_organism("Anidulans")
    for tf_gene in tf_genes_ani:
        tf_ortho_list_pox = sqlite_get_orthologs_by_src_organism("Poxalicum" , tf_gene.locus_tag)
        for tf_ortho in tf_ortho_list_pox:
            if(tf_ortho.src_gene.is_tf == 0):
                sqlite_update_gene_tf(tf_ortho.src_gene.locus_tag , 1)
        tf_ortho_list_puc = sqlite_get_orthologs_by_src_organism("Pucsensis" , tf_gene.locus_tag)
        for tf_ortho in tf_ortho_list_puc:
            if(tf_ortho.src_gene.is_tf == 0):
                sqlite_update_gene_tf(tf_ortho.src_gene.locus_tag , 1)
    # Ncrassa tf orthologs
    tf_genes_ncr = sqlite_get_tfs_by_organism("Ncrassa")
    for tf_gene in tf_genes_ncr:
        tf_ortho_list_pox = sqlite_get_orthologs_by_src_organism("Poxalicum" , tf_gene.locus_tag)
        for tf_ortho in tf_ortho_list_pox:
            if(tf_ortho.src_gene.is_tf == 0):
                sqlite_update_gene_tf(tf_ortho.src_gene.locus_tag , 1)
        tf_ortho_list_puc = sqlite_get_orthologs_by_src_organism("Pucsensis" , tf_gene.locus_tag)
        for tf_ortho in tf_ortho_list_puc:
            if(tf_ortho.src_gene.is_tf == 0):
                sqlite_update_gene_tf(tf_ortho.src_gene.locus_tag , 1)
    # Scerevisiae tf orthologs
    tf_genes_sce = sqlite_get_tfs_by_organism("Scerevisiae")
    for tf_gene in tf_genes_sce:
        tf_ortho_list_pox = sqlite_get_orthologs_by_src_organism("Poxalicum" , tf_gene.locus_tag)
        for tf_ortho in tf_ortho_list_pox:
            if(tf_ortho.src_gene.is_tf == 0):
                sqlite_update_gene_tf(tf_ortho.src_gene.locus_tag , 1)
        tf_ortho_list_puc = sqlite_get_orthologs_by_src_organism("Pucsensis" , tf_gene.locus_tag)
        for tf_ortho in tf_ortho_list_puc:
            if(tf_ortho.src_gene.is_tf == 0):
                sqlite_update_gene_tf(tf_ortho.src_gene.locus_tag , 1)

"""
    RSAT Web Service: Scan sequences for a given matrix
"""
def call_matrix_scan(service, fasta_content_str, matrix_str, uth_val, format, src_organism):

	# Wrap all arguments into a named list
	#http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
	arguments = {
                'sequence' : fasta_content_str,
                'matrix' : matrix_str,
                'matrix_format' : format,
                'uth' : ['pval '+str(uth_val)],
                'quick' : 1,
                'str' : 2,
                'origin' : 'start',
                'background_input' : 1, # this option requires 'markov'
                'markov' : 1,
                'pseudo' : 1,
                'n_treatment' : 'score',
                'background_pseudo': 0.01,
                'organism': str(src_organism)
                #'return_fields': 'pval',
 		#'verbosity' : 1
	}

	## Perform SOAP request on RSAT server
	result = service.matrix_scan(arguments)
	return result.client

"""
    Perform TFBS Predictions for each Regulatory Interaction of an organism and save in txt files
"""
def perform_tfbs_predictions(src_organism):
    tf_list = sqlite_get_tfs_by_organism(src_organism)
    pos = 0
    while (pos<len(tf_list)):
        tf_record = tf_list[pos]
        tf_locus_tag = tf_record.locus_tag
        interaction_list = sqlite_get_interactions_by_tf_locus_tag(tf_locus_tag)
        for interaction_record in interaction_list:
            tg_locus_tag = interaction_record.tg_ortho.src_gene.locus_tag
            ortho_tf_locus_tag = interaction_record.tf_ortho.ortho_gene.locus_tag
            pwm_list = sqlite_get_pwm_by_tf_locus_tag(ortho_tf_locus_tag)
            for pwm in pwm_list:
                prediction_file = out_folder + "tfbs_predictions/" + src_organism + "/" +tf_locus_tag+"-"+tg_locus_tag+"-"+pwm.motif_id+ ".txt"
                # this condition verify if exists a TFBS prediction already performed for this pwm
                # to avoid duplicated predictions with the same tg_promoter and pwm
                # note that pwms in Cis-Bp could be inferred by similarity tf_status = I
                if not os.path.isfile(prediction_file):
                    #Run RSAT
                    in_fasta_file = open(in_folder +"promoters/"+ src_organism + "/" + tg_locus_tag + "_promoter.fasta", "r")
                    fasta_content = in_fasta_file.read()
                    
                    ## Prepare matrix
                    in_matrix = in_folder + "pwms/" + pwm.motif_id + ".txt"
                    # Input matrix (file content)
                    in_matrix_file = open(in_matrix, "r")
                    matrix = in_matrix_file.read()

                    #close in files
                    in_fasta_file.close()
                    in_matrix_file.close()

                    #########################################################'
                    ## Perform SOAP request on RSAT server with the matrix and the FASTA sequence
                    lib.log.info("Call RSAT Web Service: "+'\n'+
                             " tf_locus_tag: "+ tf_locus_tag+'\n'+
                             " tg_locus_tag: "+ tg_locus_tag+'\n'+
                             " pwm_id: "+ pwm.motif_id)

                    # Call web service matrix_scan
                    if matrix == ("Pos	A	C	G	T"+'\n'):
                        lib.log.info("Transfac Matrix")
                        result = "#seq_id	ft_type	ft_name	strand	start	end	sequence	weight	Pval	ln_Pval	sig"+'\n'
                    else:
                        result = call_matrix_scan(rsat_service, fasta_content, matrix, uth_pval, matrix_format, src_organism)

                    # Output file for this matrix, from arguments in the command line
                    prediction_file = out_folder + "tfbs_predictions/" + src_organism + "/" +tf_locus_tag+"-"+tg_locus_tag+"-"+pwm.motif_id+ ".txt"

                    # Write result in output file
                    with open(prediction_file, 'w') as out_file:
                        out_file.write(result)
                    out_file.close()
                else:
                    lib.log.info("TFBS Prediction File already exists: " +tf_locus_tag+"-"+tg_locus_tag+"-"+pwm.motif_id+ ".txt")
        pos = pos + 1

"""
    Compile TFBS Predictions for each Regulatory Interaction and insert into the database
"""
def compile_tfbs_prediction(src_organism):
    lib.log.info("Compile TFBS Predictions started")
    tf_list = sqlite_get_tfs_by_organism(src_organism)
    for tf_record in tf_list:
        tf_locus_tag = tf_record.locus_tag
        interaction_list = sqlite_get_interactions_by_tf_locus_tag(tf_locus_tag)
        for interaction_record in interaction_list:
            tg_locus_tag = interaction_record.tg_ortho.src_gene.locus_tag
            ortho_tf_locus_tag = interaction_record.tf_ortho.ortho_gene.locus_tag
            pwm_list = sqlite_get_pwm_by_tf_locus_tag(ortho_tf_locus_tag)
            for pwm in pwm_list:
                prediction_file = out_folder + "tfbs_predictions/" + src_organism + "/" +tf_locus_tag+"-"+tg_locus_tag+"-"+pwm.motif_id+ ".txt"

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
                        tfbs = TFBSPrediction(pwm.motif_id, strand, start, end, sequence, weight, pval, ln_pval, sig)
                        interaction_record.add_tfbs_prediction(tfbs)
                        sqlite_insert_tfbs_prediction(interaction_record, tfbs)
                    #Create or Update Network Node
                    node = sqlite_get_network_node_by_tf_tg(tf_locus_tag, tg_locus_tag)
                    if node is None:
                        sqlite_insert_network_node(interaction_record)
                    else:
                        sqlite_update_network_node(node,interaction_record)
                in_file.close()
    lib.log.info("TFBS Predictions compilation completed")


"""
    #################### OUTPUT FILES ######################
"""    
"""
    Write TFBS Predictions Output File of an organism
"""
def write_tfbs_predictions_output_file(organism):
    lib.log.info("Writing "+organism+" TFBS Predictions Output File")
    if organism == "Poxalicum": 
        out_file_tfbs = out_file_poxalicum_tfbs
    elif organism == "Pucsensis": 
        out_file_tfbs = out_file_pucsensis_tfbs
    
    try:
        lib.log.info("Looking for TFBS Predictions")
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM tfbs_prediction WHERE src_organism = ?",[str(organism)])
        with open(out_file_tfbs, "w") as csv_file:
            csv_writer = csv.writer(csv_file, delimiter="\t")
            csv_writer.writerow([i[0] for i in cursor.description])
            csv_writer.writerows(cursor)
        csv_file.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to select data into sqlite table", error)
    lib.log.info(organism+" TFBS Predictions Output File succesfully created")

"""
    Write network output file of an organism
"""
def write_network_output_file(organism):
    lib.log.info("Writing "+organism+" Network Output File")
    
    if organism == "Poxalicum": out_file_name = out_file_poxalicum_network
    elif organism == "Pucsensis":out_file_name = out_file_pucsensis_network
    
    try:
        lib.log.info("Looking for Network Nodes.")
        cursor = sqliteConnection.cursor()
        cursor.execute("SELECT * FROM network_node WHERE organism = ?",[str(organism)])
        
        with open(out_file_name, "w") as csv_file:
            csv_writer = csv.writer(csv_file, delimiter="\t")
            csv_writer.writerow([i[0] for i in cursor.description])
            csv_writer.writerows(cursor)
        csv_file.close()
    except sqlite3.Error as error:
        lib.log.info("Failed to select data into sqlite table", error)
    lib.log.info(organism+"Network Output File successfully created")

""" 
    #################### MAIN ######################
"""
"""
    Main function of this program
""" 
if __name__ == '__main__':
    #Delete current database and create an empty new database
    #create_sqlite_database(sqlite_db_file,sqlite_db_script_file)
    
    #Create a database connection
    sqliteConnection = create_sqlite_connection(sqlite_db_file)
    
    #Execute Load Input Files
    #load_input_files()
    
    #Execute Update of TFs by orthology
    #update_tfs()
    
    #Execute Poxalicum TFBS Prediction
    #perform_tfbs_predictions("Poxalicum")
    #Execute Pucsensis TFBS Prediction
    #perform_tfbs_predictions("Pucsensis")
    
    #Compile Poxalicum TFBS Prediction
    #compile_tfbs_prediction("Poxalicum")
    #Execute Pucsensis TFBS Prediction
    #compile_tfbs_prediction("Pucsensis")
    
    #Execute Create Interactions src Pucsensis x ortho Poxalicum
    #sqlite_update_network_node_puc_by_pox_ortho()
    
    #write pucsensis output files
    write_tfbs_predictions_output_file("Poxalicum")
    write_network_output_file("Poxalicum")
    #write pucsensis output files
    write_tfbs_predictions_output_file("Pucsensis")
    write_network_output_file("Pucsensis")
    
    if (sqliteConnection):
        sqliteConnection.close()
        lib.log.info("The SQLite connection is closed")