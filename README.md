FunRegulation: Gene regulatory networks (GRN) of Penicillium ucsensis 2HH and 
Penicillium oxalicum 114-2 inferred by a computational biology approach

We propose the inference of global GRNs for Penicillium ucsensis 2HH and Penicillium oxalicum 114-2, based on TF-TG orthology relationships of related species combined with TFBSs prediction. First, global GRNs of related species (A. nidulans, N. crassa and S. cerevisiae) afford the mapping of orthologous interactions. Further, the TFBSs prediction provides accuracy to TF-TG relationships.

Python packages required:
i)   pip3 install Biopython
ii)  pip3 install suds.jurko
iii) pip3 install natsort
iv)  pip3 install pysqlrte3

  INPUT DATA:

a) Genomic data: Four fungal genomes and proteomes were used:
   - Penicillium ucsensis 2HH (this study);
   - Penicillium oxalicum 114-2 (GCA000346795.1 pdev1.0 from GenBank);
   - Aspergillus nidulansFGSC A4 (v. s10-m04-r06 from AspGD / GCF000149205.2 ASM14920v2 from GenBank);
   - Neurospora crassaOR74A (v. 12 from Broad88Institute / GCF000182925.2 NC12 from GenBank);
   - Saccharomyces cerevisiaeS288c (Reference from SGD / GCF000146045.2 R64 from GenBank).

b) Regulatory interactions from related species:
   Previous regulatory interactions of A. nidulans and N. crassa available in (Hu et al. 2018) 
   and regulatory interactions of S. cerevisiae available in YEASTRACT (Monteiro et al. 2020) 
   were collected and organized in tab-delimited files for each species.
   Regulatory interaction files:
   i)   in_file_regulatory_Ani
   ii)  in_file_regulatory_Ncr
   iii) in_file_regulatory_Sce (6 partitioned files due the data volume)

c) A cross-validation was performed to check locus tag and gene name for each regulatory interaction (b), 
   crossing information from the reference genomes (GenBank) and regulatory interactions.
   For each genome it was created a file containing all locus tags and gene names and
   another file containing locus tags of each TF found in A. nidulans, N. crassa and S. cerevisiae.
   Gene names files:
   i)   in_file_names_Puc
   ii)  in_file_names_Pox
   iii) in_file_names_Ani
   iv)  in_file_names_Ncr
   v)   in_file_names_Sce
   tf_locus_tag files:
   vi)  in_file_tfs_Ani
   vii) in_file_tfs_Ncr
   viii)in_file_tfs_Sce

d) ProteinOrtho (V6.0.15) Lechner et al. (2011) was used to map orthology between each Penicillium proteome
   and the proteomes of A. nidulans, N. crassa and S. cerevisiae, individually.
   ProteinOrtho orthology files:
   i)   in_file_ortho_Pox_Ani
   ii)  in_file_ortho_Pox_Ncr
   iii) in_file_ortho_Pox_Sce
   iv)  in_file_ortho_Puc_Pox
   v)   in_file_ortho_Puc_Ani
   v)   in_file_ortho_Puc_Ncr
   v)   in_file_ortho_Puc_Sce

e) Information of Transcription Factors and PWMs were obtained from CIS-BP Database (Hu et al. 2013).
   CIS-BP information files:
   i)   in_file_cisbp_Pox
   ii)  in_file_cisbp_Ani
   iii) in_file_cisbp_Ncr
   iv)  in_file_cisbp_Sce
  
f) PWMs in 'cis-bp' format from P.oxalicum, A. nidulans, N. crassa and S. cerevisiae 
   were obtained from CIS-BP Database (Hu et al. 2013).
   PWM files:
   i) in_folder/pwms/
  
g) Regulatory sequences of P. ucsensis 2HH and P. oxalicum 114-2 were obtained by 
   funregulation_promoter_extract.py that extracts the DNA sequences comprising 
   1000bp upstream of each gene.
   Promoter regions in fasta format:
   i)  in_folder/promoters/Poxalicum/
   ii) in_folder/promoters/Pucsensis/
  
  SQLITE3 DATABASE:

   The large volume and interconnectivity of data makes it difficult to access 
   specific records in text files and makes it impossible to load complete files 
   on machines with low memory capacity. Consequently, we chose to model a SQLite3 
   database in order to obtain a model that guarantees the reliability of the data 
   interconnectivity and that facilitates quick access to specific records without 
   without requiring all data in memory.
   SQLite3 information files:
   i) sqlite_db_folder
   ii) sqlite_db_file_name
   iii) sqlite_db_script_name
         
  OUTPUT FILES:

   Tab-delimited final output files:
   i)   out_file_pucsensis_network
   ii)  out_file_pucsensis_tfbs
   iii) out_folder/tfbs_predictions/Pucsensis/
   iv)  out_file_poxalicum_network
   i)   out_file_poxalicum_tfbs
   vi)  out_folder/tfbs_predictions/Poxalicum/


Regulatory sequences (promoter regions) of P. ucsensis 2HH and P. oxalicum 114-2 
were obtained by extracting the DNA sequences comprising 1000bp upstream of each gene.
This script uses annotation in gff3 format and whole genome sequences.

Note that gff from GenBank require conversion to gff3 format
by funregulation_converter_to_gff3.py

  INPUT FILES:

   Penicillium ucsensis 2HH (this study) genomic data:
   i)  in_file_genome (fasta format)
   ii) in_file_annotation (gff3 format)
   Penicillium oxalicum 114-2 (GCA000346795.1 pdev1.0 from GenBank) genomic data:
   i)  in_file_genome (fasta format)
   ii) in_file_annotation (gff3 format)
   
  OUTPUT FILES:
   
   Tab-delimited output file:
   i) out/extract_results/
   Promoter sequences fasta files:
   i) out/extract_misc/

Promoter sequences extraction requires annotation in GFF3 format.
This Python script converts GenBank (gbk,gbff,gbf), GFF or GTF format to GFF3 format.

  INPUT FILE:

   i) in_file_annotation (GenBank (gbk,gbff,gbf), GFF or GTF format)
   
  OUTPUT FILE:
   
   i) out_file_annotation (gff3 format)
