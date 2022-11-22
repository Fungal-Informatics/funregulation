from asyncio.subprocess import STDOUT
from subprocess import Popen, PIPE
#from django.conf import settings 
from celery import shared_task
from funRegulation.task_utils import FunRegulationBaseTask

@shared_task(bind=True, name='export_project', base=FunRegulationBaseTask)
def run_rsat():
    # RUN THE MATRIX-SCAN SCRIPT FROM RSAT
    '''
    matrix-scan -m matrixfile [-i inputfile] [-o outputfile] [-v] [-bgfile backgroundfile|-bgorder #]
        -m = Matrix file?
        -i = Arquivo de sequência
        -o = Arquivo de saída
        -v = ?
        -bgfile = This option allows to enter the background model from a background
                  model file. Background model files are tab-delimited files containing
                  the specification of oligonucleotide frequencies.
        -bgorder = ?
        command = $RSAT/perl-scripts/matrix-scan -v 1   -matrix_format transfac -m $RSAT/public_html/tmp/apache/2022/11/21/matrix-scan_2022-11-21.032537_H2IELw.matrix -pseudo 1 -decimals 1 -2str -origin end -bginput -markov 1 -bg_pseudo 0.01 -return limits -return sites -return pval -return rank -lth score 1  -uth pval 1e-4  -i $RSAT/public_html/tmp/apache/2022/11/21/tmp_sequence_2022-11-21.032537_8SojX9.fasta -seq_format fasta -n score
    '''
    rsat_path = "/usr/local/bin/matrix-scan"
    rsat_command = "-v 1 -matrix_format cis-bp -m pwm_matrix.txt -uht pval 1e-4 -i promoter.fasta -seq_format fasta -n score"
    command = ["perl", rsat_path, rsat_command]
    proc = Popen(command, stdout=PIPE, stderr=PIPE)
    output, error = proc.communicate()

    ret = proc.returncode
    if ret != 0:
        print(output)
        print('\n')
        print(error)

def run_proteinortho():
    protein_path = "/usr/local/bin/proteinortho-master/proteinortho6.pl"
    organism1 = "/home/gabriel/Downloads/TCC I/Software/ProteinOrtho/proteinortho-master/test/C.faa"
    organism2 = "/home/gabriel/Downloads/TCC I/Software/ProteinOrtho/proteinortho-master/test/E.faa"
    #
    command = ["perl", protein_path,organism1,organism2]
    
    proc = Popen(command, stdout=PIPE, stderr=PIPE)
    output, error = proc.communicate()

    ret = proc.returncode
    if ret != 0:
        print("Proteinortho failed %d %s %s" % (proc.returncode, output, error))
    else:
        print(proc.returncode, output, error)
    
    # for line in output.decode().split('\n'):
    #     print(line)

run_rsat()
#run_proteinortho()
