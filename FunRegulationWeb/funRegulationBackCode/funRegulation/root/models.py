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

class Gene:
    def __init__(self, organism, locus_tag, symbol_gene, description, is_tf):
        self.organism = organism
        self.locus_tag = locus_tag
        self.symbol_gene = symbol_gene
        self.description = description
        self.is_tf = is_tf
        
class Promoter:
    def __init__(self, locus_tag, strand, source, start, stop):
        self.locus_tag = locus_tag
        self.strand = strand
        self.source = source
        self.start = start
        self.stop = stop

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

class Orthology:
    def __init__(self, model_protein, target_protein):
        self.model_protein = model_protein
        self.target_protein = target_protein
        
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

class RegulatoryInteraction:
    def __init__(self, id, tf_locus_tag, tg_locus_tag, regulatory_function, pubmedid_source):
        self.id = id
        self.tf_locus_tag = tf_locus_tag
        self.tg_locus_tag = tg_locus_tag
        self.regulatory_function = regulatory_function
        self.pubmedid_source = pubmedid_source

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