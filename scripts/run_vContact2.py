#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    from pathlib import Path
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 

def run_vcontact2(all_vRhyme_faa, pro2viral_gn_map, tax_classification_db_dir, cluster_one_jar, outdir, threads):

    ref_viral_faa = f'{tax_classification_db_dir}/IMGVR_high-quality_phage_vOTU_representatives.faa'
    ref_pro2viral_gn_map = f'{tax_classification_db_dir}/IMGVR_high-quality_phage_vOTU_representatives_pro2viral_gn_map.csv'
    
    # Make tmp input files
    dir_path = Path(all_vRhyme_faa).parent
    
    os.system(f'cat {all_vRhyme_faa} {ref_viral_faa} > {dir_path}/combined_viral_faa.faa')
    os.system(f'cat {pro2viral_gn_map} {ref_pro2viral_gn_map} > {dir_path}/combined_pro2viral_gn_map.csv')

    # Run vcontact
    cmd = f'vcontact2 --raw-proteins {dir_path}/combined_viral_faa.faa --rel-mode Diamond --proteins-fp {dir_path}/combined_pro2viral_gn_map.csv --db None --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {cluster_one_jar} --output-dir {outdir} -t {int(threads)} -v 1> /dev/null' 
    os.system(cmd)   

all_vRhyme_faa, pro2viral_gn_map, tax_classification_db_dir, cluster_one_jar, outdir, threads = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
run_vcontact2(all_vRhyme_faa, pro2viral_gn_map, tax_classification_db_dir, cluster_one_jar, outdir, threads)  



