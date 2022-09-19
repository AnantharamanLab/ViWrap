#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    warnings.filterwarnings("ignore")   
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 

    
def run_iphop(all_vRhyme_fasta_Nlinked, iphop_outdir, iphop_db_dir, threads):
    run_cmd = f'iphop predict --fa_file {all_vRhyme_fasta_Nlinked} --out_dir {iphop_outdir} -t {threads} --db_dir {iphop_db_dir} --no_qc 1> /dev/null'
    os.mkdir(iphop_outdir)
    os.system(run_cmd)    

    
all_vRhyme_fasta_Nlinked, iphop_outdir, iphop_db_dir, threads = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
run_iphop(all_vRhyme_fasta_Nlinked, iphop_outdir, iphop_db_dir, threads)    