#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 
    
def run_vibrant(metagenomic_scaffold, viwrap_outdir, threads, virome, input_length_limit, db_dir):
    cmd = ""
    
    if virome:
        cmd = f'VIBRANT_run.py -i {metagenomic_scaffold} -folder {viwrap_outdir} -t {int(threads)} -virome -l {input_length_limit} -d {db_dir}/VIBRANT_db/databases -m {db_dir}/VIBRANT_db/files 1> /dev/null'
    else:
        cmd = f'VIBRANT_run.py -i {metagenomic_scaffold} -folder {viwrap_outdir} -t {int(threads)} -l {input_length_limit} -d {db_dir}/VIBRANT_db/databases -m {db_dir}/VIBRANT_db/files 1> /dev/null'
    os.system(cmd)     
    
metagenomic_scaffold, viwrap_outdir, threads, virome, input_length_limit, db_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]   
run_vibrant(metagenomic_scaffold, viwrap_outdir, threads, virome, input_length_limit, db_dir)    