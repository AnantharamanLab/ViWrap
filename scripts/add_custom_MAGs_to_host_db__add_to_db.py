#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings("ignore")   
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 

def add_custom_MAGs_to_host_db__add_to_db(viwrap_outdir, custom_MAGs_dir, iphop_db_dir, iphop_db_custom_dir):
    
    if custom_MAGs_dir[-1] == '/':
        custom_MAGs_dir = custom_MAGs_dir[:-1]
    
    os.system(f'iphop add_to_db --fna_dir {custom_MAGs_dir} --gtdb_dir {viwrap_outdir}/07_iPHoP_outdir/custom_MAGs_GTDB-tk_results --out_dir {iphop_db_custom_dir} --db_dir {iphop_db_dir}')
        
viwrap_outdir, custom_MAGs_dir, iphop_db_dir, iphop_db_custom_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
add_custom_MAGs_to_host_db__add_to_db(viwrap_outdir, custom_MAGs_dir, iphop_db_dir, iphop_db_custom_dir)        