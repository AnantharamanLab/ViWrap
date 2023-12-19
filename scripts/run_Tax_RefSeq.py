#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    from pathlib import Path
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call  
    warnings.filterwarnings("ignore")
    from collections import defaultdict
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 
    
    
def find_best_hits(input_diamond_outfile):
    pro2best_hit = {}

    with open(input_diamond_outfile, "r") as file:
        for line in file:
            pro, hit, bit_score = line.rstrip("\n").split("\t")
            bit_score = float(bit_score)

            if pro not in pro2best_hit or bit_score >= pro2best_hit[pro][1]:
                pro2best_hit[pro] = [hit, bit_score]

    return {pro: hit for pro, (hit, _) in pro2best_hit.items()}
   
def run_diamond_to_RefSeq_viral_protein_db(viwrap_outdir, vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, NCBI_RefSeq_viral_protein_db_dir, pro2viral_gn_map, threads, output):
    tmp_outdir = f'{viwrap_outdir}/tmp_dir_refseq'
    os.mkdir(tmp_outdir)
    
    # Step 1 Run diamond
    bin2addr = {} # bin => addr to the bin; i.e., vRhyme_bin_10 => path/to/the/dir/vRhyme_bin_10.faa
    walk = os.walk(vRhyme_best_bin_dir)
    for path, dir_list, file_list in walk:
        for file_name in file_list:
            if "faa" in file_name: 
                file_name_with_path = os.path.join(path, file_name)
                bin_name = Path(file_name_with_path).stem
                bin2addr[bin_name] = file_name_with_path
                
    walk2 = os.walk(vRhyme_unbinned_viral_gn_dir)
    for path, dir_list, file_list in walk2:
        for file_name in file_list:
            if "faa" in file_name: 
                file_name_with_path = os.path.join(path, file_name)
                bin_name = Path(file_name_with_path).stem
                bin2addr[bin_name] = file_name_with_path  

    diamond_cmd = []
    for bin_name in bin2addr:
        bin_addr = bin2addr[bin_name]
        each_cmd = f'diamond blastp -q {bin_addr} -p 1 --db {NCBI_RefSeq_viral_protein_db_dir}/NCBI_RefSeq_viral.dmnd --evalue 0.00001 --query-cover 50 --subject-cover 50 -k 10000 -o {tmp_outdir}/{bin_name}.diamond_out.txt -f 6 --quiet 1> /dev/null'
        diamond_cmd.append(each_cmd)
    
    n = int(threads) # The number of parallel processes
    for j in range(max(int(len(diamond_cmd)/n + 1), 1)):
        procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in diamond_cmd[j*n: min((j+1)*n, len(diamond_cmd))] ]
        for p in procs:
            p.wait()     

    # Step 2 Summarize the result            
    # Store 2.1 Store pro information in a bin
    pro2bin = {} # pro (header wo arrow) => bin_name
    bin2pro_num = {} # bin_name => pro_num; Store the number of proteins in each bin
    
    with open(pro2viral_gn_map, "r") as lines:
        for line in lines:
            line = line.rstrip("\n")
            if not line.startswith("protein_id"):
                pro, bin_name = line.split(",")[0], line.split(",")[1]
                pro2bin[pro] = bin_name
                bin2pro_num[bin_name] = bin2pro_num.get(bin_name, 0) + 1
    lines.close()
    
    # Store 2.2 Store the diamond db pro 2 tax info
    NCBI_RefSeq_viral_protein2tax = {} # pro => tax
    with open(f'{NCBI_RefSeq_viral_protein_db_dir}/pro2ictv_8_rank_tax.txt' ,"r") as lines:
        for line in lines:
            line = line.rstrip("\n")
            pro, tax = line.split("\t", 1)[0], line.split("\t", 1)[1]
            NCBI_RefSeq_viral_protein2tax[pro] = tax

    # Store 2.3 Store the best hits and to see whether >= 30% of the proteins for a bin have a hit to Viral RefSeq
    bin2best_hits = defaultdict(list) # bin_name => [best_hits]
    # Only record this if >= 30% of the proteins for a faa have a hit to Viral RefSeq
    for bin_name in bin2addr:
        diamond_out = f'{tmp_outdir}/{bin_name}.diamond_out.txt'
        if os.path.exists(diamond_out):
            pro2best_hit = find_best_hits(f'{tmp_outdir}/{bin_name}.diamond_out.txt')
            bin_involved = set() # Store the bin that have sequences inside have the best hit; now the bin is the old name
            for pro in pro2best_hit:
                bin_name3 = pro2bin[pro] # the bin_name3 is the old name
                bin_involved.add(bin_name3)
                
            # Split the best hit result into each bin
            for bin_name2 in bin_involved: # This bin_name2 is the subset, only for bin_involved
                pro2best_hit_in_this_bin = {} # # Store all the proteins that have best hits in this bin
                    
                for pro in pro2best_hit:
                    if pro2bin[pro] == bin_name2:
                        pro2best_hit_in_this_bin[pro] = pro2best_hit[pro]
                            
                pro_num_w_best_hit = len(pro2best_hit_in_this_bin) # The number of proteins within in bin have the best hits
                bin_pro_num = bin2pro_num[bin_name2] # The total protein name from this bin
                if float(pro_num_w_best_hit/bin_pro_num) >= 0.3: # To see if >=30% of the proteins for a bin have a hit to Viral RefSeq
                    for pro in pro2best_hit_in_this_bin:
                        best_hit = pro2best_hit_in_this_bin[pro]
                        bin2best_hits[bin_name2].append(best_hit)
                
    # Store 2.4 Get the consensus affiliation based on the best hits of individual proteins (>= 50 majority rule)
    bin2consensus_tax = {} # bin => consensus_tax
    for bin_name in bin2best_hits:
        best_hits = bin2best_hits[bin_name]
        
        tax2freq = {} # The frequency of tax; tax => frequency
        for best_hit in best_hits:
            tax = NCBI_RefSeq_viral_protein2tax[best_hit]
            tax2freq[tax] = tax2freq.get(tax, 0) + 1
                
        tax_w_highest_freq = "" 
        highest_freq = 0
        for tax in tax2freq:
            freq = tax2freq[tax]
            if freq > highest_freq:
                tax_w_highest_freq = tax
                highest_freq = freq
                
        prec_of_tax_w_highest_freq = 0
        prec_of_tax_w_highest_freq = highest_freq / len(best_hits)
        
        if float(prec_of_tax_w_highest_freq) >= 0.5: # Only store the consensus tax if there is
            bin2consensus_tax[bin_name] = tax_w_highest_freq
           
    # Step 2.5 Remove tmp dir 
    os.system(f'rm -r {tmp_outdir}')
    
    # Step 2.6 Write to output
    f = open(output, 'w')
    for bin_name in bin2consensus_tax:
        f.write(f'{bin_name}\t{bin2consensus_tax[bin_name]}\n')
    f.close()    
    
viwrap_outdir, vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, NCBI_RefSeq_viral_protein_db_dir, pro2viral_gn_map, threads, output = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]
run_diamond_to_RefSeq_viral_protein_db(viwrap_outdir, vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, NCBI_RefSeq_viral_protein_db_dir, pro2viral_gn_map, threads, output)    