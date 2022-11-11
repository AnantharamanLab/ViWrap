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
    
def run_virsorter2_1st(metagenomic_scaffold, virsorter_outdir, threads, input_length_limit):
    vs_cmd = f"virsorter run --keep-original-seq -i {metagenomic_scaffold} -w {virsorter_outdir}/pass1 --min-length {input_length_limit} --min-score 0.5 -j {threads} all 1> /dev/null" 
    os.system(vs_cmd)    
      
metagenomic_scaffold, virsorter_outdir, threads, input_length_limit = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] 
run_virsorter2_1st(metagenomic_scaffold, virsorter_outdir, threads, input_length_limit)

  