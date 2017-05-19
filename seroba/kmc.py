import os
from seroba import common, external_progs

class Error (Exception): pass

ext_progs = external_progs.ExternalProgs()


def run_kmc(sequence_file,kmer_size,temp_dir,prefix):
    file_format = common.detect_sequence_format(sequence_file)
    kmer_count_files = os.path.join(temp_dir,prefix)
    if file_format == 'fa':
        param = ' -ci1 -m1 -t1 -fm'
    elif file_format =='fq':
        param = ' -ci4  -m1 -t1'
    print(param)
    command1=[ext_progs.exe('kmc') + ' -k'+kmer_size,param,sequence_file,kmer_count_files,temp_dir]
    print(' '.join(command1))
    os.system(' '.join(command1))
    return kmer_count_files

def run_kmc_intersect(kmer_count_files,temp_dir,db):
    temp_inter=os.path.join(temp_dir,'inter')
    command2 = [ext_progs.exe('kmc_tools'), 'simple',kmer_count_files,db,'intersect',temp_inter]
    print(' '.join(command2))
    os.system(' '.join(command2))

    return temp_inter

def run_kmc_hist(temp_inter,temp_dir):
    temp_hist=os.path.join(temp_dir,'hist')
    command3 = [ext_progs.exe('kmc_tools'),'transform',temp_inter,'histogram',temp_hist]
    os.system(' '.join(command3))
    return temp_hist
