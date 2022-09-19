#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    from pathlib import Path
    from glob import glob
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call      
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 
    
def get_hmmsearch_result(hmmsearch_result):
    pro2info = {} # pro => [query, score, evalue]
    with open(hmmsearch_result, 'r') as lines:
        for line in lines:
            if not line.startswith('#'):
                line = line.rstrip("\n")
                line = re.sub(' +', ' ', line)
                tmp = line.split(' ')
                pro, query, score, evalue = tmp[0], tmp[2], tmp[5], tmp[4]
                pro2info[pro] = [query, score, evalue] 
    lines.close()            
    return pro2info         
 
def run_hmmsearch_to_marker_VOG_HMM_db(vog_marker_table, viwrap_outdir, vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, tax_classification_db_dir, pro2viral_gn_map, threads, output):
    tmp_outdir = f'{viwrap_outdir}/tmp_dir_vog'
    os.mkdir(tmp_outdir)
    
    # Step 1 Run hmmsearch
    bin2addr = {} # bin => addr to the bin; i.e., vRhyme_bin_10 => path/to/the/dir/vRhyme_bin_10.faa
    file_names1 = glob(f'{vRhyme_best_bin_dir}/*.faa')
    file_names2 = glob(f'{vRhyme_unbinned_viral_gn_dir}/*.faa')
    file_names = file_names1 + file_names2
    for file_name in file_names:
        bin_name = Path(file_name).stem
        bin2addr[bin_name] = file_name
  
    hmmsearch_cmd = []
    for bin_name in bin2addr:
        bin_addr = bin2addr[bin_name]
        each_cmd = f'hmmsearch -E 0.01 --cpu 1 --tblout {tmp_outdir}/{bin_name}.hmmsearch_result.txt {tax_classification_db_dir}/marker_VOG.hmm {bin_addr} 1> /dev/null'
        hmmsearch_cmd.append(each_cmd)
    
    n = int(threads) # The number of parallel processes
    for j in range(max(int(len(hmmsearch_cmd)/n + 1), 1)):
        procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in hmmsearch_cmd[j*n: min((j+1)*n, len(hmmsearch_cmd))] ]
        for p in procs:
            p.wait()        

    # Step 2 Get marker VOG info
    vog_marker_list = {} # vog => tax
    with open(vog_marker_table, "r") as lines:
        for line in lines:
            line = line.rstrip("\n")
            if not line.startswith("#"):
                vog, tax = line.split("\t")[0], line.split("\t")[2]
                vog_marker_list[vog] = tax
    lines.close()

    # Step 3 Filter hmmsearch result to get protein hits to VOG marker hash
    pro2vog = {}
    hmmsearch_results = glob(f'{tmp_outdir}/*.hmmsearch_result.txt')
    for hmmsearch_result in hmmsearch_results: 
        pro2info = get_hmmsearch_result(hmmsearch_result)
        for pro in pro2info: 
            vog, score, evalue = pro2info[pro][0], pro2info[pro][1], pro2info[pro][2]
            if float(score) >= 40 and float(evalue) <= 0.00001:
                pro2vog[pro] = vog
     
    # Step 4 Find a consensus taxonomy for each bin (simple plurality rule)
    pro2bin = {} # pro (header wo arrow) => bin_name
    with open(pro2viral_gn_map, "r") as lines:
        for line in lines:
            line = line.rstrip("\n")
            if not line.startswith("protein_id"):
                pro, bin_name = line.split(",")[0], line.split(",")[1]
                pro2bin[pro] = bin_name
    lines.close()           
           
    bin2pro_hits = {} # bin_name => [pro_hits]; Store bin with pro hits
    for pro in pro2vog:
        bin_name = pro2bin[pro]
        if bin_name not in bin2pro_hits:
            bin2pro_hits[bin_name] = [pro]
        else:
            bin2pro_hits[bin_name].append(pro)
            
    bin2consensus_tax = {} # bin_name => consensus_tax (simple plurality rule)
    for bin_name in bin2pro_hits:
        pro_hits = bin2pro_hits[bin_name]
        
        tax2freq = {}
        for pro in pro_hits:
            vog = pro2vog[pro]
            tax = vog_marker_list[vog]
            tax2freq[tax] = tax2freq.get(tax, 0) + 1
                
        consensus_tax = ""
        consensus_tax_freq = 0
        for tax in tax2freq:
            if tax2freq[tax] > consensus_tax_freq:
                consensus_tax_freq = tax2freq[tax]
                consensus_tax = tax
        
        bin2consensus_tax[bin_name] = consensus_tax

    # Step 5 Remove tmp dir 
    os.system(f'rm -r {tmp_outdir}')
     
    # Step 6 Write to output
    f = open(output, 'w')
    for bin_name in bin2consensus_tax:
        f.write(f'{bin_name}\t{bin2consensus_tax[bin_name]}\n')
    f.close()

    
vog_marker_table, viwrap_outdir, vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, tax_classification_db_dir, pro2viral_gn_map, threads, output = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8]
run_hmmsearch_to_marker_VOG_HMM_db(vog_marker_table, viwrap_outdir, vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, tax_classification_db_dir, pro2viral_gn_map, threads, output)    