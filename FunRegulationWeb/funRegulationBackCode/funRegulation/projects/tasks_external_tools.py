from asyncio.subprocess import STDOUT
from subprocess import Popen, PIPE
import urllib.parse
from root import general_functions as LenzSoftware
#from django.conf import settings 
#from celery import shared_task
#from funRegulation.task_utils import FunRegulationBaseTask

#@shared_task(bind=True, name='export_project', base=FunRegulationBaseTask)
def run_rsat():
    rsat_pwms = "/home/gabriel/packages/rsat/perl-scripts/M01890_2.00.txt"
    rsat_promoter = "/home/gabriel/packages/rsat/perl-scripts/FOXG_10713_promoter.fasta"
    
    cmd = "matrix-scan -v 1 -matrix_format cis-bp -m "+ rsat_pwms+" -pseudo 1 -decimals 1 -2str -origin end -bginput -markov 1 -bg_pseudo 0.01 -return limits -return sites -return pval -return rank -lth score 1 -uth pval 1e-4 -i "+rsat_promoter+" -seq_format fasta -n score"
    
    rsat_call = Popen(cmd, shell=True,stdout=PIPE,stderr=PIPE)
    out, error = rsat_call.communicate()
    out_folder = '/home/gabriel/Desktop/FunRegulationBack-end/funregulation/FunRegulationWeb/funRegulationBackCode/funRegulation/projects/out_files/'
    out_file = out_folder + "temp.txt"

    ret = rsat_call.returncode
    if ret != 0:
        print("RSAT failed %d %s %s" % (rsat_call.returncode, out, error))
    else:
        LenzSoftware.construct_grn_tfbs_predictions('5',out)
        # with open(out_file,'w') as teste:
        #     teste.write(str(out))
        # teste.close()
        # with open (out_file) as in_file:
        #     for line in in_file:
        #         if line.startswith("#"): 
        #             continue
        #         parts = line.strip().split("\\t")
        #         strand = urllib.parse.unquote(parts[76])
        #         start = urllib.parse.unquote(parts[77])
        #         end = urllib.parse.unquote(parts[78])
        #         sequence = urllib.parse.unquote(parts[79])
        #         weight = urllib.parse.unquote(parts[80])
        #         pval = urllib.parse.unquote(parts[81])
        #         ln_pval = urllib.parse.unquote(parts[82])
        #         sig = urllib.parse.unquote(parts[83])
        #     in_file.close()
        # print(strand, start, end, sequence, weight,pval,ln_pval,sig)

def run_proteinortho():
    protein_path = "/usr/local/bin/proteinortho-master/proteinortho6.pl"
    organism1 = "/home/gabriel/Downloads/TCC I/Software/ProteinOrtho/proteinortho-master/test/C.faa"
    organism2 = "/home/gabriel/Downloads/TCC I/Software/ProteinOrtho/proteinortho-master/test/E.faa"
    
    command = ["perl", protein_path,organism1,organism2]
    
    proc = Popen(command, stdout=PIPE, stderr=PIPE)
    output, error = proc.communicate()

    ret = proc.returncode
    if ret != 0:
        print("Proteinortho failed %d %s %s" % (proc.returncode, output, error))
    else:
        for line in output.decode().split('\n'):
            print(line)
        #print(proc.returncode, output, error)
    
    # for line in output.decode().split('\n'):
    #     print(line)

run_rsat()
#run_proteinortho()
