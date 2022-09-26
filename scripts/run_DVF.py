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
    
def run_dvf(metagenomic_scaffold, dvf_outdir, input_length_limit, db_dir):
    cmd = f"dvf.py -i {metagenomic_scaffold} -o {dvf_outdir} -l {int(input_length_limit)} -m {db_dir} 1> /dev/null"
    os.system(cmd) 
    
metagenomic_scaffold, dvf_outdir, input_length_limit, db_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]   
run_dvf(metagenomic_scaffold, dvf_outdir, input_length_limit, db_dir)    

