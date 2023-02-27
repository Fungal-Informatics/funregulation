import psycopg2
import lib.library as lib
import os, platform
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import urllib.parse
from root.models import *
from collections import namedtuple

#Initialize Folder Paths
in_folder = ""
out_folder = ""

#Initialize file Paths
in_file_organisms = os.path.join(in_folder,'')
in_file_genes = os.path.join(in_folder,'')
in_file_proteins = os.path.join(in_folder,'')
in_file_pwms = os.path.join(in_folder,'')
in_file_model_regulatory = os.path.join(in_folder,'')
in_file_orthology = os.path.join(in_folder,'')
in_file_genome = os.path.join(in_folder,'')

upstream = -1000
downstream = 0

dbConnection = None

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

def parse_protein_file(in_file_proteins):
    lib.log.info("Parsing "+ in_file_proteins)
    for rec in SeqIO.parse(in_file_proteins, 'fasta'):
        
        #when locus_tag != protein_id
        #rec.description = re.search(r'gene:(.*?) transcript:', rec.description).group(1)
        #protein = Protein(rec.description, rec.id,'','','','','','','','','','')
        
        #when locus_tag == protein_id
        protein = Protein(rec.id, rec.id,'','','','','','','','','','')
        insert_protein(protein)

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
        #lib.log.info(source)
        lib.log.info(protein_id)

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

def construct_grn_tfbs_predictions(organism_id, rsat_result):
    tf_list = list()
    pwm_list = list()
    regulatory_interactions_list = list()
    
    # parse fasta file and turn into dictionary
    genome = SeqIO.to_dict(SeqIO.parse(open(in_file_genome), 'fasta'))
    
    tf_list = select_tfs_by_organism(organism_id)
    for tf in tf_list:
        pwm_list = select_pwms_by_locus_tag(tf.locus_tag)
        for pwm in pwm_list:
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
                        #result = call_matrix_scan(rsat_service, promoter_sequence, matrix)
                        result = rsat_result
                        # Write result in output file
                        with open(prediction_file, 'w') as out_file:
                            out_file.write(result)
                        out_file.close()

                        with open (prediction_file) as in_file:
                            for line in in_file:
                                if line.startswith("#"): 
                                    continue
                                parts = line.strip().split("\\t")
                                #Create new TFBS prediction for each RSAT prediction result
                                strand = urllib.parse.unquote(parts[76])
                                start = urllib.parse.unquote(parts[77])
                                end = urllib.parse.unquote(parts[78])
                                sequence = urllib.parse.unquote(parts[79])
                                weight = urllib.parse.unquote(parts[80])
                                pval = urllib.parse.unquote(parts[81])
                                ln_pval = urllib.parse.unquote(parts[82])
                                sig = urllib.parse.unquote(parts[83])
                                tfbs = TFBS(0, regulatory_interaction.id, pwm.id, strand, start, end, sequence, weight, pval, ln_pval, sig)
                                insert_tfbs_prediction(tfbs)
                            in_file.close()
                else:
                    lib.log.info("TFBS Prediction File already exists: " +regulatory_interaction.tf_locus_tag+"-"+regulatory_interaction.tg_locus_tag+"-"+pwm.motif_id+ ".txt")
