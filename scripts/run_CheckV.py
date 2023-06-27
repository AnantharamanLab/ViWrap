#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    warnings.filterwarnings("ignore")
    from pathlib import Path
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call    
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 
    
def run_checkv(input_dir, outdir, threads, checkv_db_dir):
    checkv_cmd = []
    walk = os.walk(input_dir)
    for path, dir_list, file_list in walk:
        for file_name in file_list:
            if "fasta" in file_name:
                file_name_with_path = os.path.join(path, file_name)
                file_name_stem = Path(file_name).stem
                each_cmd_1 = f'checkv completeness {file_name_with_path} {outdir}/{file_name_stem} -t 1 -d {checkv_db_dir} 1> /dev/null'
                each_cmd_2 = f'checkv complete_genomes {file_name_with_path} {outdir}/{file_name_stem} -d {checkv_db_dir} 1> /dev/null'
                each_cmd_3 = f'checkv quality_summary {file_name_with_path} {outdir}/{file_name_stem} -d {checkv_db_dir} 1> /dev/null'
                each_cmd = each_cmd_1 + ';' + each_cmd_2 + ';' + each_cmd_3
                checkv_cmd.append(each_cmd)
                
    n = int(threads) # The number of parallel processes
    for j in range(max(int(len(checkv_cmd)/n + 1), 1)):
        procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in checkv_cmd[j*n: min((j+1)*n, len(checkv_cmd))] ]
        for p in procs:
            p.wait()       
    
input_dir, outdir, threads, checkv_db_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
run_checkv(input_dir, outdir, threads, checkv_db_dir)    