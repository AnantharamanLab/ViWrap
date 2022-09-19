#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from glob import glob
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call  
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 
    
def run_drep(dRep_outdir, viral_genus_genome_list_dir, threads, dRep_length_limit):
    dRep_cmd = []
    viral_genus_genome_lists = glob(f'{viral_genus_genome_list_dir}/viral_genus_genome_list.*.txt')
    viral_genus_genome_lists_non_singleton = []
    for viral_genus_genome_list in viral_genus_genome_lists:
        with open(viral_genus_genome_list, 'r') as fp:
            line_num = len(fp.readlines())
            if line_num != 1:
                viral_genus_genome_lists_non_singleton.append(viral_genus_genome_list)

    for viral_genus_genome_list in viral_genus_genome_lists_non_singleton:
        VC = Path(viral_genus_genome_list).stem.split(".")[1]
        each_cmd = f'dRep dereplicate {dRep_outdir}/Output.{VC} -p 1 -g {viral_genus_genome_list} -l {dRep_length_limit} --ignoreGenomeQuality -pa 0.8 -sa 0.95 -nc 0.85 -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -centW 0 1> /dev/null'
        dRep_cmd.append(each_cmd)

    n = int(threads) # The number of parallel processes
    for j in range(max(int(len(dRep_cmd)/n + 1), 1)):
        procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in dRep_cmd[j*n: min((j+1)*n, len(dRep_cmd))] ]
        for p in procs:
            p.wait()      
    
dRep_outdir, viral_genus_genome_list_dir, threads, dRep_length_limit = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
run_drep(dRep_outdir, viral_genus_genome_list_dir, threads, dRep_length_limit)    