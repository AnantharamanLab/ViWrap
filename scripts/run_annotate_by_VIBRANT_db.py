#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    from pathlib import Path  
    import math
    from collections import defaultdict  
    from glob import glob    
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call  
    from Bio import SeqIO    
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 
    
def store_seq(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict
    
    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.rstrip("\n") # Remove "\n" in the end
            if ">" in line:
                if " " in line or "\t" in line: # Break at the first " " or "\t"
                    spliter = ""
                    for i in range(len(line)):
                        if line[i] == " " or line[i] == "\t":
                            spliter = line[i]
                            break 
                           
                    head = line.split(f'{spliter}', 1)[0]
                    seq_dict[head] = ""
                else:
                    head = line
                    seq_dict[head] = ""
            else:
                seq_dict[head] += line
            
    seq_lines.close()
    
    return seq_dict

def write_down_seq(seq_dict, path_to_file): 
    # Two inputs are required:
    # (1) The dict of the sequence
    # (2) The path that you want to write your sequence down
    
    seq_file = open(path_to_file,"w")
    for head in seq_dict:
        seq_file.write(head + "\n")
        seq_file.write(seq_dict[head] + "\n")
    seq_file.close()   

def chuncker(list_to_split, chunk_size):
    list_of_chunks =[]
    start_chunk = 0
    end_chunk = start_chunk+chunk_size
    while end_chunk <= len(list_to_split)+chunk_size:
        chunk_ls = list_to_split[start_chunk: end_chunk]
        list_of_chunks.append(chunk_ls)
        start_chunk = start_chunk +chunk_size
        end_chunk = end_chunk+chunk_size    
    return list_of_chunks    
    
def split_seq(input_seq, split_num, output_seq_folder):
    # Step 1 Store the seq dict
    input_seq_dict = store_seq(input_seq)
    
    # Step 2 Make list of each seq dict
    input_seq_dict_keys_list = list(input_seq_dict.keys())
    chunk_size = math.ceil(len(input_seq_dict_keys_list) / int(split_num))
    input_seq_dict_keys_list_of_chunks = chuncker(input_seq_dict_keys_list, chunk_size)
    
    list_of_seq_dicts = [] # [seq_dicts]
    for chunk in input_seq_dict_keys_list_of_chunks:
        seq_dict = {}
        for header in chunk:
            seq_dict[header] = input_seq_dict[header]
        list_of_seq_dicts.append(seq_dict)

    # Step 3 Write down individual seq files
    stem_name = Path(input_seq).stem
    suffix = Path(input_seq).suffix
    if os.path.exists(output_seq_folder):
        sys.exit(f"The output dir - {output_seq_folder} - for storing splited fasta/faa files has been created!")
    else:
        os.mkdir(output_seq_folder)
    
    for i in range(len(list_of_seq_dicts)):
        j = i + 1
        seq_dict = list_of_seq_dicts[i]
        output_seq_file = os.path.join(output_seq_folder, f"{stem_name}.chunk_{j}{suffix}")
        if len(seq_dict) > 0:
            write_down_seq(seq_dict, output_seq_file)
        
def get_hmmsearch_result(hmmsearch_result):
    pro2info = {} # pro => [query, query_accession, evalue, score]
    with open(hmmsearch_result, 'r') as lines:
        for line in lines:
            if not line.startswith('#'):
                line = line.rstrip("\n")
                line = re.sub(' +', ' ', line)
                tmp = line.split(' ')
                pro, query, query_accession, evalue, score = tmp[0], tmp[2], tmp[3], tmp[4], tmp[5]
                pro2info[pro] = [query, query_accession, evalue, score]
    lines.close()            
    return pro2info   

## The function for obtaining geNomad ffn file 
def extract_genes(fasta_file, faa_file, output_file):
    # Parse the gene coordinates, strand information, and full headers from the faa file
    gene_info = {}
    with open(faa_file, 'r') as faa:
        for record in SeqIO.parse(faa, "fasta"):
            parts = record.description.split('#')
            if len(parts) >= 4:
                start = int(parts[1].strip())
                end = int(parts[2].strip())
                strand = int(parts[3].strip())
                full_header = record.description
                gene_id = record.id
                gene_info[gene_id] = (start, end, strand)

    # Extract gene sequences from the fasta file using the coordinates and strand information
    extracted_genes = []
    with open(fasta_file, 'r') as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            for gene_id, (start, end, strand) in gene_info.items():
                # Compare the modified gene_id with the record.id
                if gene_id.rsplit('_', 1)[0] == record.id:
                    if '|provirus_' in record.id: # Dealing with provirus coordinates - the coordinates are from the original scaffold
                        provirus_start = int(record.id.split(' ')[0].split('|provirus_')[1].split('_', 1)[0])
                        start_real = start - provirus_start + 1
                        end_real = end - provirus_start + 1
                        if strand == 1:
                            gene_sequence = record.seq[start_real-1:end_real]
                        else:
                            # For reverse strand, reverse complement the sequence
                            gene_sequence = record.seq[start_real-1:end_real].reverse_complement()                        
                    else:
                        if strand == 1:
                            gene_sequence = record.seq[start-1:end]
                        else:
                            # For reverse strand, reverse complement the sequence
                            gene_sequence = record.seq[start-1:end].reverse_complement()
                    extracted_genes.append((gene_id, gene_sequence))

    # Write the extracted gene sequences to the output file with the original headers
    with open(output_file, 'w') as output:
        for gene_id, gene in extracted_genes:
            output.write(f">{gene_id}\n{gene}\n")
    
def run_annotate_by_vibrant_db(VIBRANT_db, identify_method, virsorter_outdir, dvf_outdir, genomad_outdir, out_dir, threads):
    final_virus_fasta_file = ''
    KEGG_hmm_file = os.path.join(VIBRANT_db, 'databases/KEGG_profiles_prokaryotes.HMM')
    Pfam_hmm_file = os.path.join(VIBRANT_db, 'databases/Pfam-A_v32.HMM')
    VOG_hmm_file = os.path.join(VIBRANT_db, 'databases/VOGDB94_phage.HMM')
    
    if identify_method == 'vs':
        final_virus_fasta_file = os.path.join(virsorter_outdir, 'final_vs2_virus.fasta')
    elif identify_method == 'dvf':
        final_virus_fasta_file = os.path.join(dvf_outdir, 'final_dvf_virus.fasta')
    elif identify_method == 'genomad':
        final_virus_fasta_file = os.path.join(genomad_outdir, 'final_genomad_virus.fasta')        

    # Step 1 Prodigal annotate fasta file to get faa and ffn files or extract ffn from fasta by faa
    final_virus_faa_file = final_virus_fasta_file.replace('.fasta', '.faa', 1)
    final_virus_ffn_file = final_virus_fasta_file.replace('.fasta', '.ffn', 1)
    if identify_method == 'vs' or identify_method == 'dvf':
        prodigal_cmd = f"prodigal -i {final_virus_fasta_file} -a {final_virus_faa_file} -d {final_virus_ffn_file} -p meta -q -o /dev/null"
        os.system(prodigal_cmd) 
    elif identify_method == 'genomad':     
        extract_genes(final_virus_fasta_file, final_virus_faa_file, final_virus_ffn_file)        
        
    # Step 2 Get all split faa addresses
    output_seq_folder = os.path.join(out_dir, 'tmp_dir_split_faa')
    split_seq(final_virus_faa_file, threads, output_seq_folder)    
    all_faa_addrs = glob(f"{output_seq_folder}/*.faa")

    # Step 3 Run hmmsearch against KEGG database
    tmp_dir_kegg_hmmsearch_results = os.path.join(out_dir, 'tmp_dir_kegg_hmmsearch_results')
    if os.path.exists(tmp_dir_kegg_hmmsearch_results):
        sys.exit(f"The output dir - {tmp_dir_kegg_hmmsearch_results}  - for storing KEGG hmmseach results has been created!")
    else:
        os.mkdir(tmp_dir_kegg_hmmsearch_results)
    
    kegg_hmmsearch_cmds = []
    for faa_addr in all_faa_addrs:
        faa_stem = Path(faa_addr).stem
        kegg_hmmtbl = os.path.join(tmp_dir_kegg_hmmsearch_results, f"{faa_stem}.KEGG.hmmtbl")
        kegg_temp = os.path.join(tmp_dir_kegg_hmmsearch_results, f"{faa_stem}_temp.txt")
        each_cmd = f"hmmsearch --tblout {kegg_hmmtbl} --noali -T 40 --cpu 1 -o {kegg_temp} {KEGG_hmm_file} {faa_addr}"
        kegg_hmmsearch_cmds.append(each_cmd)
    
    n = int(threads) # The number of parallel processes
    for j in range(max(int(len(kegg_hmmsearch_cmds)/n + 1), 1)):
        procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in kegg_hmmsearch_cmds[j*n: min((j+1)*n, len(kegg_hmmsearch_cmds))] ]
        for p in procs:
            p.wait()  

    # Step 4 Run hmmsearch against Pfam database
    tmp_dir_pfam_hmmsearch_results = os.path.join(out_dir, 'tmp_dir_pfam_hmmsearch_results')
    if os.path.exists(tmp_dir_pfam_hmmsearch_results):
        sys.exit(f"The output dir - {tmp_dir_pfam_hmmsearch_results}  - for storing Pfam hmmseach results has been created!")
    else:
        os.mkdir(tmp_dir_pfam_hmmsearch_results)
    
    pfam_hmmsearch_cmds = []
    for faa_addr in all_faa_addrs:
        faa_stem = Path(faa_addr).stem
        pfam_hmmtbl = os.path.join(tmp_dir_pfam_hmmsearch_results, f"{faa_stem}.Pfam.hmmtbl")
        pfam_temp = os.path.join(tmp_dir_pfam_hmmsearch_results, f"{faa_stem}_temp.txt")
        each_cmd = f"hmmsearch --tblout {pfam_hmmtbl} --noali -T 40 --cpu 1 -o {pfam_temp} {Pfam_hmm_file} {faa_addr}"
        pfam_hmmsearch_cmds.append(each_cmd)
    
    n = int(threads) # The number of parallel processes
    for j in range(max(int(len(pfam_hmmsearch_cmds)/n + 1), 1)):
        procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in pfam_hmmsearch_cmds[j*n: min((j+1)*n, len(pfam_hmmsearch_cmds))] ]
        for p in procs:
            p.wait()  

    # Step 5 Run hmmsearch against VOG database
    tmp_dir_vog_hmmsearch_results = os.path.join(out_dir, 'tmp_dir_vog_hmmsearch_results')
    if os.path.exists(tmp_dir_vog_hmmsearch_results):
        sys.exit(f"The output dir - {tmp_dir_vog_hmmsearch_results}  - for storing VOG hmmseach results has been created!")
    else:
        os.mkdir(tmp_dir_vog_hmmsearch_results)
    
    vog_hmmsearch_cmds = []
    for faa_addr in all_faa_addrs:
        faa_stem = Path(faa_addr).stem
        vog_hmmtbl = os.path.join(tmp_dir_vog_hmmsearch_results, f"{faa_stem}.VOG.hmmtbl")
        vog_temp = os.path.join(tmp_dir_vog_hmmsearch_results, f"{faa_stem}_temp.txt")
        each_cmd = f"hmmsearch --tblout {vog_hmmtbl} --noali -T 40 --cpu 1 -o {vog_temp} {VOG_hmm_file} {faa_addr}"
        vog_hmmsearch_cmds.append(each_cmd)
    
    n = int(threads) # The number of parallel processes
    for j in range(max(int(len(vog_hmmsearch_cmds)/n + 1), 1)):
        procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in vog_hmmsearch_cmds[j*n: min((j+1)*n, len(vog_hmmsearch_cmds))] ]
        for p in procs:
            p.wait()  
            
    # Step 6 Parse hmmsearch results
        #KEGG -> query
        #Pfam -> query_accession
        #VOG -> query  
    ## Step 6.1 Parse KEGG hmmsearch results
    KEGG_hmm_result = {} # pro => [query, evalue, score]
    KEGG_hmmtbls = glob(os.path.join(tmp_dir_kegg_hmmsearch_results, '*.KEGG.hmmtbl'))
    for hmmtbl in KEGG_hmmtbls:
        pro2info = get_hmmsearch_result(hmmtbl)
        for pro in pro2info:
            query, evalue, score = pro2info[pro][0], pro2info[pro][2], pro2info[pro][3]
            KEGG_hmm_result[pro] = [query, evalue, score]
            
    ## Step 6.2 Parse Pfam hmmsearch results
    Pfam_hmm_result = {} # pro => [query_accession, evalue, score]
    Pfam_hmmtbls = glob(os.path.join(tmp_dir_pfam_hmmsearch_results, '*.Pfam.hmmtbl'))
    for hmmtbl in Pfam_hmmtbls:
        pro2info = get_hmmsearch_result(hmmtbl)
        for pro in pro2info:
            query_accession, evalue, score = pro2info[pro][1], pro2info[pro][2], pro2info[pro][3]
            Pfam_hmm_result[pro] = [query_accession, evalue, score]

    ## Step 6.3 Parse VOG hmmsearch results
    VOG_hmm_result = {} # pro => [query, evalue, score]
    VOG_hmmtbls = glob(os.path.join(tmp_dir_vog_hmmsearch_results, '*.VOG.hmmtbl'))
    for hmmtbl in VOG_hmmtbls:
        pro2info = get_hmmsearch_result(hmmtbl)
        for pro in pro2info:
            query, evalue, score = pro2info[pro][0], pro2info[pro][2], pro2info[pro][3]
            VOG_hmm_result[pro] = [query, evalue, score]            
            
    ## Step 6.4 Store KO, Pfam, VOG info
    KO2info = {} # KO => [AMG, KO name, v_score]
    Pfam2info = {} # Pfam => [Pfam name, v_score]
    VOG2info = {} # VOG => [VOG name, v_score]
    
    AMG_KO = [] # [KO] Store the AMG KOs
    with open(os.path.join(VIBRANT_db, 'files/VIBRANT_AMGs.tsv'),'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('KO'):
                KO = line
                AMG_KO.append(KO)
    lines.close()
    
    ## Store "VIBRANT_categories.tsv" to get vscore for KO, Pfam, and VOG
    hmm2v_score = {} # hmm => v_score (hmm is either KO, Pfam, or VOG)
    with open(os.path.join(VIBRANT_db, 'files/VIBRANT_categories.tsv'),'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('id'):
                tmp = line.split('\t')
                hmm = tmp[0]
                v_score = float(tmp[1]) / 100
                hmm2v_score[hmm] = str(v_score)

    with open(os.path.join(VIBRANT_db, 'files/VIBRANT_names.tsv'),'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            hmm, name = tmp[0], tmp[1]
            if hmm.startswith('VOG'):
                v_score = ''
                if hmm in hmm2v_score:
                    v_score = hmm2v_score[hmm]
                VOG2info[hmm] = [name, v_score]                
            elif hmm.startswith('K'):
                v_score = ''
                if hmm in hmm2v_score:
                    v_score = hmm2v_score[hmm]                                
                if hmm in AMG_KO:
                    KO2info[hmm] = ['AMG', name, v_score]
                else:
                    KO2info[hmm] = ['', name, v_score]
            else:
                v_score = ''
                if hmm in hmm2v_score:
                    v_score = hmm2v_score[hmm]                    
                Pfam2info[hmm] = [name, v_score]     
    lines.close()            
                
    # Step 7 Write down annotation result
    annotation_file = ''
    if identify_method == 'vs':
        annotation_file = os.path.join(virsorter_outdir, 'final_vs2_virus.annotation.txt')
    elif identify_method == 'dvf':
        annotation_file = os.path.join(dvf_outdir, 'final_dvf_virus.annotation.txt')
    elif identify_method == 'genomad':
        annotation_file = os.path.join(genomad_outdir, 'final_genomad_virus.annotation.txt')
        
    f = open(annotation_file, 'w')
    header_list = ['protein', 'scaffold', 'KO', 'AMG', 'KO name', 'KO evalue', 'KO score', 'KO v-score', 'Pfam', 'Pfam name', 'Pfam evalue', 'Pfam score', 'Pfam v-score', 'VOG', 'VOG name', 'VOG evalue', 'VOG score', 'VOG v-score']    
    header = '\t'.join(header_list)
    f.write(header + '\n')
    all_pro_seq = {}
    for faa_addr in all_faa_addrs:
        pro_seq = store_seq(faa_addr)
        all_pro_seq.update(pro_seq)
    
    for pro_w_array in all_pro_seq:
        pro = pro_w_array.replace('>', '' , 1)
        scf = pro.rsplit('_', 1)[0]
        KO, AMG, KO_name, KO_evalue, KO_score, KO_v_score = '', '', '', '', '', ''
        if pro in KEGG_hmm_result:
            KO, KO_evalue, KO_score = KEGG_hmm_result[pro][0], KEGG_hmm_result[pro][1], KEGG_hmm_result[pro][2]
            AMG, KO_name, KO_v_score = KO2info[KO][0], KO2info[KO][1], KO2info[KO][2]
        Pfam, Pfam_name, Pfam_evalue, Pfam_score, Pfam_v_score = '', '', '', '', ''
        if pro in Pfam_hmm_result:
            Pfam, Pfam_evalue, Pfam_score = Pfam_hmm_result[pro][0], Pfam_hmm_result[pro][1], Pfam_hmm_result[pro][2]
            Pfam_name, Pfam_v_score = Pfam2info[Pfam][0], Pfam2info[Pfam][1]
        VOG, VOG_name, VOG_evalue, VOG_score, VOG_v_score = '', '', '', '', ''
        if pro in VOG_hmm_result:
            VOG, VOG_evalue, VOG_score = VOG_hmm_result[pro][0], VOG_hmm_result[pro][1], VOG_hmm_result[pro][2]
            VOG_name, VOG_v_score = VOG2info[VOG][0], VOG2info[VOG][1]
        line_list =  [pro, scf, KO, AMG, KO_name, KO_evalue, KO_score, KO_v_score, Pfam, Pfam_name, Pfam_evalue, Pfam_score, Pfam_v_score, VOG, VOG_name, VOG_evalue, VOG_score, VOG_v_score]
        line = '\t'.join(line_list) 
        f.write(line + '\n')  
    f.close()   

    # Step 7 Remove all tmp dirs    
    os.system(f"rm -rf {output_seq_folder} {tmp_dir_kegg_hmmsearch_results} {tmp_dir_pfam_hmmsearch_results} {tmp_dir_vog_hmmsearch_results}")
               
    
VIBRANT_db, identify_method, virsorter_outdir, dvf_outdir, genomad_outdir, out_dir, threads = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]
run_annotate_by_vibrant_db(VIBRANT_db, identify_method, virsorter_outdir, dvf_outdir, genomad_outdir, out_dir, threads)       