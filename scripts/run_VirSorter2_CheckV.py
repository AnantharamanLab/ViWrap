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
    
def run_virsorter2_checkv(virsorter_outdir, threads, checkv_db_dir):

    checkv_cmd = f"checkv end_to_end {virsorter_outdir}/pass/final-viral-combined.fa {virsorter_outdir}/CheckV_result -t {threads} -d {checkv_db_dir} 1> /dev/null"
    os.system(checkv_cmd)   

    cat_cmd = f"cat {virsorter_outdir}/CheckV_result/proviruses.fna {virsorter_outdir}/CheckV_result/viruses.fna > {virsorter_outdir}/CheckV_result/combined.fna"
    os.system(cat_cmd) 
     
virsorter_outdir, threads, checkv_db_dir = sys.argv[1], sys.argv[2], sys.argv[3]
run_virsorter2_checkv(virsorter_outdir, threads, checkv_db_dir)

  