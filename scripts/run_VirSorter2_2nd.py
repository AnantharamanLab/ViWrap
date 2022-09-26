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
    
def run_virsorter2_2nd(virsorter_outdir, threads, input_length_limit):
    vs_2nd_cmd = f"virsorter run --seqname-suffix-off --viral-gene-enrich-off --prep-for-dramv -i {virsorter_outdir}/CheckV_result_1st/combined.fna -w {virsorter_outdir}/pass2 --min-length {input_length_limit} --min-score 0.5 -j {threads} all 1> /dev/null"
    os.system(vs_2nd_cmd)    

virsorter_outdir, threads, input_length_limit = sys.argv[1], sys.argv[2], sys.argv[3]
run_virsorter2_2nd(virsorter_outdir, threads, input_length_limit)

  