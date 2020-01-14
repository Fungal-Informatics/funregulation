#!/usr/bin/env python

#Wrapper script for Funannotate package.
import pkg_resources
import sys
import os
import subprocess
import inspect

version = 0.1 #pkg_resources.require("funregulation")[0].version
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

def flatten(l):
    flatList = []
    for elem in l:
        # if an element of a list is a list
        # iterate over this list and add elements to flatList 
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList

default_help = """
Usage:       funregulation <command> <arguments>
version:     %s
Description: Funregulation is a pipeline to extract intergenic regions and annotate promoters and transcription binding sites.
    
Command:     promoterextract     Run promoter extraction pipeline
             #tbspredict         Run transcription binding sites prediction pipeline
             promoterannotate    Assign annotation to promoter elements
             #tbsannotate        Assign annotation to transcription binding sites
             
Written by Alexandre Lenz (2017-2021) arlenz@ucs.br / alenz@uneb.br
        """ % version

if len(sys.argv) > 1:
    
    if sys.argv[1] == 'promoterextract':
        help = """
Usage:       funregulation %s <arguments>
version:     %s
Description: Script takes genome multi-fasta file and GFF3 annotation file to do a whole
             genome intergenic regions extraction using the start_codon. 
    
Required:  -g, --genome           Genome multi-fasta file.
           -a, --annotation       GFF3 annotation file.
           -o, --out              Output folder name.
            
Written by Alexandre Lenz (2017-2019) arlenz@ucs.br / alenz@uneb.br
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funregulation-promoterextract.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print (help)
            sys.exit(1)
    elif sys.argv[1] == 'promoterannotate':
        help = """
Usage:       funregulation %s <arguments>
version:     %s
Description: Script annotates the elements from funregulation promoterextract.  It pulls
             annotation info using  
             ElemeNT: A Computational Tool for Detecting Core Promoter Elements
             http://www.ncbi.nlm.nih.gov/pubmed/26226151

Required:    -i, --input        Folder from funregulation promoter extract
          
Written by Alexandre Lenz (2017-2019) arlenz@ucs.br / alenz@uneb.br
        """ % (sys.argv[1], version)
        
    elif sys.argv[1] == 'compare':
        help = """
Usage:       funregulation %s <arguments>
version:     %s
Description: Script does comparative regulatory elements between funregulation genomes.  Output
             is graphs, CSV files --> visualized in web-browser.  
    
Required:    -i, --input         List of funregulation genome folders
Optional:    -o, --out           Output folder name. Default: funregulation_compare
             
Written by Alexandre Lenz (2017-2019) arlenz@ucs.br / alenz@uneb.br
        """ % (sys.argv[1], version)
       
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funregulation-compare.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print (help)
            sys.exit(1)
else:
    print (default_help)