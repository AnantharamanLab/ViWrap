#!/usr/bin/env python3

''' 
Aim: Store frequently used python3 functions
'''

try:
    import warnings
    import sys
    import os
    import re
    import pandas as pd
    from statistics import mean
    from collections import defaultdict
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from glob import glob
    import pyfastx # For fastq and fasta reading and parsing
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
                if (" " or "\t") in line: # Break at the first " " or "\t"
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
    
def store_seq_with_full_head(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict
    
    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.replace("\n","") # Remove "\n"
            if line[0] == ">":
                head = line
                seq_dict[head] = ""
            else:
                seq_dict[head] += line
            
    seq_lines.close()
    
    return seq_dict
    
def get_gene_seq(input_gene_file): # Add the filename to the header; 
    # The input sequence file should be a file with full path
    head = "" # Store the header line
    gene_dict = {} # Store the sequence dict
    filename = Path(input_gene_file).stem # Get the stem name from a file with full path
    
    with open(input_gene_file,"r") as gene_lines:
        for line in gene_lines:
            line = line.replace("\n","") # Remove "\n"
            if ">" in line[0]:
                if "\s" in line:
                    head = ">" + filename + "~~" + line.replace(">", "").split("\s")[0]
                else:
                    head = ">" + filename + "~~" + line.replace(">", "")
                gene_dict[head] = ""    
            else: 
                gene_dict[head] += line
                gene_dict[head].replace("*","")
    
    gene_lines.close()
            
    return gene_dict   

def write_down_seq(seq_dict, path_to_file): 
    # Two inputs are required:
    # (1) The dict of the sequence
    # (2) The path that you want to write your sequence down
    
    seq_file = open(path_to_file,"w")
    for head in seq_dict:
        seq_file.write(head + "\n")
        seq_file.write(seq_dict[head] + "\n")
    seq_file.close()       
    
def make_unbinned_viral_gn(viral_scaffold, vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir):
    viral_scaffold_faa = viral_scaffold.rsplit(".", 1)[0] + ".faa"
    viral_scaffold_ffn = viral_scaffold.rsplit(".", 1)[0] + ".faa"       
    viral_scaffold_fasta_dict = store_seq_with_full_head(viral_scaffold)
    viral_scaffold_faa_dict = store_seq_with_full_head(viral_scaffold_faa)
    viral_scaffold_ffn_dict = store_seq_with_full_head(viral_scaffold_ffn)
    
    # Step 1 Make unnbinned fasta dict and write down unbinned viral genome fasta file
    viral_scaffold_fasta_binned_dict = {} # scaffold_id (NODE_10610_length_8667_cov_0.658730) => bin_name (vRhyme_10)

    walk = os.walk(vRhyme_best_bin_dir)
    for path,dir_list,file_list in walk:
        for file_name in file_list:
            if "fasta" in file_name: 
                file_name_with_path = os.path.join(path, file_name)
                with open(file_name_with_path, "r") as fasta_lines:
                    for line in fasta_lines:
                        if ">" in line:
                            line = line.rstrip("\n")
                            bin_name = line.replace(">", "", 1).split("__", 1)[0]
                            scaffold_id = line.replace(">", "", 1).split("__", 1)[1]
                            viral_scaffold_fasta_binned_dict[scaffold_id] = bin_name
                fasta_lines.close() 

    viral_scaffold_fasta_unbinned_dict = {} # scaffold_id (NODE_10610_length_8667_cov_0.658730) => unbinned_gn_name (vRhyme_unbinned_1)
    i = 1 
    for header in viral_scaffold_fasta_dict:
        scaffold_id = header.replace(">", "", 1)
        if scaffold_id not in viral_scaffold_fasta_binned_dict:
            viral_scaffold_fasta_unbinned_dict[scaffold_id] = f'vRhyme_unbinned_{i}'
            i += 1
    
    os.mkdir(vRhyme_unbinned_viral_gn_dir)    
    
    for scaffold_id in viral_scaffold_fasta_unbinned_dict:
        unbinned_fasta_file = open(f'{vRhyme_unbinned_viral_gn_dir}/{viral_scaffold_fasta_unbinned_dict[scaffold_id]}.fasta',"w")
        unbinned_fasta_file.write(f'>{viral_scaffold_fasta_unbinned_dict[scaffold_id]}__{scaffold_id}\n')
        scaffold_id_w_array = ">" + scaffold_id
        unbinned_fasta_file.write(f'{viral_scaffold_fasta_dict[scaffold_id_w_array]}\n')
        unbinned_fasta_file.close()
        
    # Step 2 Write down unbinned viral genome faa file    
    viral_scaffold_faa_unbinned_dict = {} # pro_id (NODE_10811_length_8534_cov_0.218494_1  (1..228)        1       PF02229.16      Transcriptional Coactivator p15 (PC4)) => unbinned_gn_name (vRhyme_unbinned_1)
    for header in viral_scaffold_faa_dict:
        pro_id = header.replace(">", "", 1)
        scaffold_id = pro_id.split("\t", 1)[0].rsplit("_", 1)[0]
        if scaffold_id in viral_scaffold_fasta_unbinned_dict:
            viral_scaffold_faa_unbinned_dict[pro_id] = viral_scaffold_fasta_unbinned_dict[scaffold_id]
   
    
    for scaffold_id in viral_scaffold_fasta_unbinned_dict:
        unbinned_faa_file = open(f'{vRhyme_unbinned_viral_gn_dir}/{viral_scaffold_fasta_unbinned_dict[scaffold_id]}.faa',"w")
        for pro_id in viral_scaffold_faa_unbinned_dict:
            scaffold_id_within = pro_id.split("\t", 1)[0].rsplit("_", 1)[0]
            if (scaffold_id_within == scaffold_id):
                unbinned_faa_file.write(f'>{viral_scaffold_faa_unbinned_dict[pro_id]}__{pro_id}\n')
                pro_id_w_array = ">" + pro_id
                unbinned_faa_file.write(f'{viral_scaffold_faa_dict[pro_id_w_array]}\n')
        unbinned_faa_file.close()    
        
    # Step 3 Write down unbinned viral genome ffn file    
    viral_scaffold_ffn_unbinned_dict = {} # pro_id (NODE_10811_length_8534_cov_0.218494_1  (1..228)        1       PF02229.16      Transcriptional Coactivator p15 (PC4)) => unbinned_gn_name (vRhyme_unbinned_1)
    for header in viral_scaffold_ffn_dict:
        pro_id = header.replace(">", "", 1)
        scaffold_id = pro_id.split("\t", 1)[0].rsplit("_", 1)[0]
        if scaffold_id in viral_scaffold_fasta_unbinned_dict:
            viral_scaffold_ffn_unbinned_dict[pro_id] = viral_scaffold_fasta_unbinned_dict[scaffold_id]
    
    for scaffold_id in viral_scaffold_fasta_unbinned_dict:
        unbinned_ffn_file = open(f'{vRhyme_unbinned_viral_gn_dir}/{viral_scaffold_fasta_unbinned_dict[scaffold_id]}.ffn',"w")
        for pro_id in viral_scaffold_ffn_unbinned_dict:
            scaffold_id_within = pro_id.split("\t", 1)[0].rsplit("_", 1)[0]
            if (scaffold_id_within == scaffold_id):
                unbinned_ffn_file.write(f'>{viral_scaffold_ffn_unbinned_dict[pro_id]}__{pro_id}\n')
                pro_id_w_array = ">" + pro_id
                unbinned_ffn_file.write(f'{viral_scaffold_ffn_dict[pro_id_w_array]}\n')
        unbinned_ffn_file.close()  
        
    # Step 4 Write down best bin faa file
    viral_scaffold_faa_binned_dict = {} # pro_id (NODE_10811_length_8534_cov_0.218494_1  (1..228)        1       PF02229.16      Transcriptional Coactivator p15 (PC4)) => unbinned_gn_name (vRhyme_bin_1)
    for header in viral_scaffold_faa_dict:
        pro_id = header.replace(">", "", 1)
        scaffold_id = pro_id.split("\t", 1)[0].rsplit("_", 1)[0]
        if scaffold_id in viral_scaffold_fasta_binned_dict:
            viral_scaffold_faa_binned_dict[pro_id] = viral_scaffold_fasta_binned_dict[scaffold_id]

    binned_faa_name2scaffolds = {} # binned_faa_name to scaffolds
    for scaffold_id in viral_scaffold_fasta_binned_dict:
        binned_faa_name = viral_scaffold_fasta_binned_dict[scaffold_id].replace("vRhyme","vRhyme_bin",1)
        if binned_faa_name not in binned_faa_name2scaffolds:
            binned_faa_name2scaffolds[binned_faa_name] = scaffold_id
        else:
            binned_faa_name2scaffolds[binned_faa_name] += "\t" + scaffold_id

    for binned_faa_name in binned_faa_name2scaffolds:
        binned_faa_file = open(f'{vRhyme_best_bin_dir}/{binned_faa_name}.faa',"w") 
        scaffolds = binned_faa_name2scaffolds[binned_faa_name].split("\t")
        for scaffold_id in scaffolds:
            for pro_id in viral_scaffold_faa_binned_dict:
                scaffold_id_within = pro_id.split("\t", 1)[0].rsplit("_", 1)[0]
                if (scaffold_id_within == scaffold_id):
                    binned_faa_file.write(f'>{viral_scaffold_faa_binned_dict[pro_id]}__{pro_id}\n')
                    pro_id_w_array = ">" + pro_id
                    binned_faa_file.write(f'{viral_scaffold_faa_dict[pro_id_w_array]}\n')
        binned_faa_file.close()
        
    # Step 5 Write down best bin ffn file
    viral_scaffold_ffn_binned_dict = {} # pro_id (NODE_10811_length_8534_cov_0.218494_1  (1..228)        1       PF02229.16      Transcriptional Coactivator p15 (PC4)) => unbinned_gn_name (vRhyme_bin_1)
    for header in viral_scaffold_ffn_dict:
        pro_id = header.replace(">", "", 1)
        scaffold_id = pro_id.split("\t", 1)[0].rsplit("_", 1)[0]
        if scaffold_id in viral_scaffold_fasta_binned_dict:
            viral_scaffold_ffn_binned_dict[pro_id] = viral_scaffold_fasta_binned_dict[scaffold_id]

    binned_ffn_name2scaffolds = {} # binned_ffn_name to scaffolds
    for scaffold_id in viral_scaffold_fasta_binned_dict:
        binned_ffn_name = viral_scaffold_fasta_binned_dict[scaffold_id].replace("vRhyme","vRhyme_bin",1)
        if binned_ffn_name not in binned_ffn_name2scaffolds:
            binned_ffn_name2scaffolds[binned_ffn_name] = scaffold_id
        else:
            binned_ffn_name2scaffolds[binned_ffn_name] += "\t" + scaffold_id

    for binned_ffn_name in binned_ffn_name2scaffolds:
        binned_ffn_file = open(f'{vRhyme_best_bin_dir}/{binned_ffn_name}.ffn',"w") 
        scaffolds = binned_ffn_name2scaffolds[binned_ffn_name].split("\t")
        for scaffold_id in scaffolds:
            for pro_id in viral_scaffold_ffn_binned_dict:
                scaffold_id_within = pro_id.split("\t", 1)[0].rsplit("_", 1)[0]
                if (scaffold_id_within == scaffold_id):
                    binned_ffn_file.write(f'>{viral_scaffold_ffn_binned_dict[pro_id]}__{pro_id}\n')
                    pro_id_w_array = ">" + pro_id
                    binned_ffn_file.write(f'{viral_scaffold_ffn_dict[pro_id_w_array]}\n')
        binned_ffn_file.close()         
    
def get_pro2viral_gn_map(vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, pro2viral_gn_map):
    pro2viral_gn_dict = {}
    walk = os.walk(vRhyme_best_bin_dir)
    for path,dir_list,file_list in walk:
        for file_name in file_list:
            if "faa" in file_name: 
                viral_gn = Path(file_name).stem
                file_name_with_path = os.path.join(path, file_name)
                with open(file_name_with_path, "r") as faa_lines:
                    for line in faa_lines:
                        if ">" in line:
                            line = line.rstrip("\n")
                            pro = line.replace(">", "", 1).split("\t", 1)[0]
                            pro2viral_gn_dict[pro] = viral_gn
                faa_lines.close()          
                            
    walk2 = os.walk(vRhyme_unbinned_viral_gn_dir)
    for path, dir_list, file_list in walk2:
        for file_name in file_list:
            if "faa" in file_name:
                viral_gn = Path(file_name).stem
                file_name_with_path = os.path.join(path, file_name)
                with open(file_name_with_path,"r") as faa_lines2:
                    for line in faa_lines2:
                        if ">" in line:
                            line = line.rstrip("\n")
                            pro = line.replace(">", "", 1).split("\t", 1)[0]
                            pro2viral_gn_dict[pro] = viral_gn
                faa_lines2.close()            
                            
    file = open(pro2viral_gn_map,"w")
    file.write('protein_id,contig_id,keywords\n')
    for pro in pro2viral_gn_dict:
        # print(pro)
        # print(pro2viral_gn_dict[pro])
        file.write(f'{pro},{pro2viral_gn_dict[pro]},None\n')
    file.close()
   
def combine_all_vRhyme_faa(vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, all_vRhyme_faa):
    walk = os.walk(vRhyme_best_bin_dir)
    walk2 = os.walk(vRhyme_unbinned_viral_gn_dir)
    all_vRhyme_faa_dict = {}
    
    for path, dir_list, file_list in walk:
        for file_name in file_list:
            if "faa" in file_name:
                file_name_with_path = os.path.join(path, file_name)
                each_vRhyme_faa_dict = store_seq(file_name_with_path)
                all_vRhyme_faa_dict.update(each_vRhyme_faa_dict) # Add dict from one to another one

    for path, dir_list, file_list in walk2:
        for file_name in file_list:
            if "faa" in file_name:
                file_name_with_path = os.path.join(path, file_name)
                each_vRhyme_faa_dict = store_seq(file_name_with_path)
                all_vRhyme_faa_dict.update(each_vRhyme_faa_dict) # Add dict from one to another one
                
    write_down_seq(all_vRhyme_faa_dict, all_vRhyme_faa)       

def combine_all_vRhyme_fasta(vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, all_vRhyme_fasta):
    walk = os.walk(vRhyme_best_bin_dir)
    walk2 = ''
    if vRhyme_unbinned_viral_gn_dir:
        walk2 = os.walk(vRhyme_unbinned_viral_gn_dir)
    all_vRhyme_fasta_dict = {}
    
    for path, dir_list, file_list in walk:
        for file_name in file_list:
            if "fasta" in file_name:
                file_name_with_path = os.path.join(path, file_name)
                each_vRhyme_faa_dict = store_seq(file_name_with_path)
                all_vRhyme_fasta_dict.update(each_vRhyme_faa_dict) # Add dict from one to another one

    if walk2:
        for path, dir_list, file_list in walk2:
            for file_name in file_list:
                if "fasta" in file_name:
                    file_name_with_path = os.path.join(path, file_name)
                    each_vRhyme_faa_dict = store_seq(file_name_with_path)
                    all_vRhyme_fasta_dict.update(each_vRhyme_faa_dict) # Add dict from one to another one
                
    write_down_seq(all_vRhyme_fasta_dict, all_vRhyme_fasta)     
   
def get_genus_cluster_info(genome_by_genome_file, genus_cluster_info):
    genus_cluster_dict = {} # VC => VC, all gn
    all_gns = set()
    clustered_gn = set()
    
    with open(genome_by_genome_file,"r") as file_lines:
        for line in file_lines:
            line = line.rstrip('\n')
            if not line.startswith('Genome,'):
                gn = line.split(",")[0]
                if 'vRhyme' in gn:
                    all_gns.add(gn)
                    genus_confidence_score = line.split(",")[9]
                    VC = line.split(",")[3]
                    if genus_confidence_score != '':
                        if float(genus_confidence_score) >= 0.9:
                            clustered_gn.add(gn)
                            if VC not in genus_cluster_dict:
                                genus_cluster_dict[VC] = [VC, gn]
                            else:
                                genus_cluster_dict[VC][1] += ";" + gn                  
    file_lines.close()

    i = 1
    for gn in all_gns:
        if gn not in clustered_gn:
            VC = f'UnclusteredGenus_{i}'
            genus_cluster_dict[VC] = [VC, gn]
            i += 1
    
    # Write down genus cluster info result
    file = open(genus_cluster_info, "w")
    file.write('#VC,genomes\n') 
    for VC in genus_cluster_dict:
        gns = genus_cluster_dict[VC][1]
        file.write(f'{VC},{gns}\n')           
    file.close()
    
def Nlinker(infolder, outdir, extension, n):
# Copied from vRhyme auxiliary scripts by Kristopher Kieft, UW-Madison
# Using n of Ns to link scaffolds
    files = os.listdir(infolder)
    len_ext = len(extension)
    files = [i for i in files if i[-len_ext:] == extension]

    if len(files) == 0:
        sys.stderr.write("\nError: No input files were identified. Verify that the input folder and extension are correct. Exiting.\n\n")
        exit()

    N_string = ''.join("N" * n)
    for f in files:
        file = infolder + "/" + f
        base = f.rsplit(".",1)[0]
        with open(file, 'r') as fasta:
            seq = ''
            for line in fasta:
                if not line.startswith(">"):
                    seq += line.strip("\n")
                    seq += N_string

        with open(outdir + "/" + base + '.linked.' + extension, "w") as outfile:
            outfile.write(">" + base + "\n" + seq[:-n] + "\n") # -n to remove last n added characters         
                
def parse_checkv_result(input_dir, outfile):
    walk = os.walk(input_dir)
    header = "" # Store the header line of quality summary file
    checkv_summary = [] # Store all checkv_summary_result 
    for path, dir_list, file_list in walk:
        for dir_name in dir_list:
            dir_name_with_path = os.path.join(path, dir_name)
            
            # Get only the first level dir
            if 'tmp' not in dir_name_with_path: 
                each_quality_summary_file = f'{dir_name_with_path}/quality_summary.tsv'
                with open(each_quality_summary_file, "r") as each_file:
                    for line in each_file:
                        line = line.rstrip("\n")
                        if line.startswith("contig_id"):
                            header = line
                        else:
                            checkv_summary.append(line)
                each_file.close()            
                        
    f = open(outfile, "w")
    f.write(header + "\n")
    for line in checkv_summary:
        f.write(line + "\n")
    f.close()    
    
def get_gn_list_for_genus(genus_cluster_info, dRep_outdir, vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir):
    genus_dict = {} # VC => gns
    with open(genus_cluster_info, "r") as genus_cluster:
        for line in genus_cluster:
            line = line.rstrip("\n")
            if 'VC' in line.split(",", 1)[0] and line[0] != '#':
                VC = line.split(",", 1)[0]
                genus_dict[VC] = line.split(",", 1)[1]
    genus_cluster.close()            

    gn_address = {} # gn stem name => the full path to each genome
    walk = os.walk(vRhyme_best_bin_dir)
    for path, dir_list, file_list in walk:
        for file_name in file_list:
            if "fasta" in file_name:
                file_name_with_path = os.path.join(path, file_name)
                file_name_stem = Path(file_name).stem
                gn_address[file_name_stem] = file_name_with_path
   
    walk2 = os.walk(vRhyme_unbinned_viral_gn_dir)
    for path, dir_list, file_list in walk2:
        for file_name in file_list:
            if "fasta" in file_name:
                file_name_with_path = os.path.join(path, file_name)
                file_name_stem = Path(file_name).stem
                gn_address[file_name_stem] = file_name_with_path

    os.mkdir(dRep_outdir)   
    os.mkdir(f'{dRep_outdir}/viral_genus_genome_list')      
    
    for VC in genus_dict:
        f = open(f'{dRep_outdir}/viral_genus_genome_list/viral_genus_genome_list.{VC}.txt', "w")
        gns = genus_dict[VC].split(";")
        for gn in gns:
            gn_w_full_path = gn_address[gn]
            f.write(gn_w_full_path + "\n")
        f.close()              
    
def parse_dRep(viwrap_outdir, dRep_outdir, species_cluster_info, genus_cluster_info, viral_genus_genome_list_dir):
    gn2VC = {} # gn => VC or UnclusteredGenus
    with open(genus_cluster_info, "r") as genus_cluster:
        for line in genus_cluster:
            line = line.rstrip("\n")
            if line[0] != '#':
                VC = line.split(",", 1)[0]
                gns = line.split(",", 1)[1].split(";")
                for gn in gns:
                    gn2VC[gn] = VC                
    genus_cluster.close()   
    
    species_cluster_dict = {} # species_rep => species_rep, gns, genus

    walk = os.walk(dRep_outdir)
    for path, dir_list, file_list in walk:   
        for dir_name in dir_list:
            if 'Output' in dir_name:
                dir_name_with_path = os.path.join(path, dir_name)
                Cdb = dir_name_with_path + "/data_tables/Cdb.csv"
                Wdb = dir_name_with_path + "/data_tables/Wdb.csv"
                Bdb = dir_name_with_path + "/data_tables/Bdb.csv"
                
                cluster2species_rep = {} # cluster => species_rep (within this genus)
                gn2cluster = {}
                
                if os.path.exists(Cdb) and os.path.exists(Wdb):
                    with open(Wdb, "r") as Wdb_file:
                        for line in Wdb_file:
                            line = line.rstrip("/n")
                            if line.split(",", 1)[0] != 'genome':
                                species_rep = line.split(",", 1)[0].split(".", 1)[0]
                                cluster = line.split(",")[1]
                                cluster2species_rep[cluster] = species_rep
                    Wdb_file.close()            
                                
                    with open(Cdb, "r") as Cdb_file:
                        for line in Cdb_file:
                            line = line.rstrip("/n")
                            if line.split(",", 1)[0] != 'genome':
                                gn = line.split(",", 1)[0].split(".", 1)[0]
                                cluster = line.split(",")[1]
                                gn2cluster[gn] = cluster
                    Cdb_file.close()  
                else:
                    with open(Bdb, "r") as Bdb_file:
                        for line in Bdb_file:
                            line = line.rstrip("/n")
                            if line.split(",", 1)[0] != 'genome':
                                species_rep = line.split(",", 1)[0].split(".", 1)[0]
                                cluster = line.split(",")[1]
                                cluster2species_rep[cluster] = species_rep
                                gn2cluster[species_rep] = cluster
                    Bdb_file.close()                      
                    
                for cluster in cluster2species_rep:
                    species_rep = cluster2species_rep[cluster]
                    species_cluster_dict[species_rep] = ['', '', '']
                    species_cluster_dict[species_rep][0] = species_rep
                    species_cluster_dict[species_rep][2] = gn2VC[species_rep]
                        
                    gns = [] # Store the genomes within this cluster
                    for gn in gn2cluster:
                        if gn2cluster[gn] == cluster:
                            gns.append(gn)
                            
                    species_cluster_dict[species_rep][1] = ';'.join(gns)       
                                
    viral_genus_genome_lists = glob(f'{viral_genus_genome_list_dir}/viral_genus_genome_list.*.txt')
    viral_genus_genome_lists_singleton = []
    for viral_genus_genome_list in viral_genus_genome_lists:
        with open(viral_genus_genome_list, 'r') as fp:
            line_num = len(fp.readlines())
            if line_num == 1:
                viral_genus_genome_lists_singleton.append(viral_genus_genome_list)
    
    for viral_genus_genome_list in viral_genus_genome_lists_singleton:
        with open(viral_genus_genome_list, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                singeton_gn = Path(line).stem.split('.')[0]
                species_cluster_dict[singeton_gn] = [singeton_gn, singeton_gn, gn2VC[singeton_gn]]              
              
    # Store UnclusteredGenus
    for gn in gn2VC:
        if 'UnclusteredGenus' in gn2VC[gn]:
            species_cluster_dict[gn] = [gn, gn, gn2VC[gn]]
            
    # Write down species_cluster_info
    f = open(f'{viwrap_outdir}/Species_cluster_info.txt',"w")
    f.write('#species_rep,genomes,genus\n')
    for species_rep in species_cluster_dict:
        genus = species_cluster_dict[species_rep][2]
        gns = species_cluster_dict[species_rep][1]
        f.write(f'{species_rep},{gns},{genus}\n')
    f.close()    

def get_virus_raw_abundance(mapping_result_dir, vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, virus_raw_abundance):
    # Step 1 Get gn2scaffolds dict
    gn2scaffolds = {} # gn => [scaffolds]
    
    all_gn_addrs1 = glob(f'{vRhyme_best_bin_dir}/*.fasta')
    all_gn_addrs2 = glob(f'{vRhyme_unbinned_viral_gn_dir}/*.fasta')
    all_gn_addrs = all_gn_addrs1 + all_gn_addrs2
    
    for gn_adds in all_gn_addrs: 
        gn = Path(gn_adds).stem
        seqs = store_seq(gn_adds)
        scaffolds = [x.replace('>', '', 1).split('__', 1)[1] for x in seqs]
        gn2scaffolds[gn] = scaffolds
        
    # Step 2 Store coverm raw coverage table  
    coverm_raw_table = pd.read_csv(f'{mapping_result_dir}/all_coverm_raw_result.txt', sep = '\t', index_col = 0)
    coverm_raw_dict = coverm_raw_table.to_dict() # col => row => value
    
    # Step 3 Get gn coverage
    ## Step 3.1 Get bams list
    bams = []
    for col in list(coverm_raw_dict.keys()):
        if '.bam' in col and '.bam-var' not in col:
            bams.append(col)
    
    ## Step 3.2 Get gn2bam2coverage dict
    gn2bam2coverage = defaultdict(dict) # gn => bam => coverage
    for gn in gn2scaffolds:
        for bam in bams:
            coverages = [] # Store all scaffold coverages
            scaffolds = gn2scaffolds[gn]
            for scaffold in scaffolds:
                if '_fragment_' in scaffold:
                    scaffold = scaffold.rsplit('_fragment_', 1)[0]
                if '||' in scaffold:
                    scaffold = scaffold.rsplit('||', 1)[0]
                coverage = coverm_raw_dict[bam][scaffold]
                coverages.append(coverage)
            gn_coverage = mean(coverages)
            gn2bam2coverage[gn][bam] = gn_coverage
            
    # Step 4 Write down dict
    gn2bam2coverage_df = pd.DataFrame(gn2bam2coverage).T
    gn2bam2coverage_df.columns = gn2bam2coverage_df.columns.str.replace('.filtered.bam', '')
    gn2bam2coverage_df.fillna(0, inplace=True)
    gn2bam2coverage_df.to_csv(virus_raw_abundance, sep='\t')

def get_read_info(metaG_reads):
    # (1) Test if the read pair is of the same size
    # (2) Get the info of read pair: read_count, read_total_base, average_length
    
    metaG_reads_list = metaG_reads.split(',')
    sample2read_info = {} # sample => [read_count, read_base]
    
    if len(metaG_reads_list) / 2 == 1:
        sample = Path(metaG_reads_list[0]).stem.rsplit('_', 1)[0]
        
        fq1 = pyfastx.Fastq(metaG_reads_list[0])
        fq2 = pyfastx.Fastq(metaG_reads_list[1])

        fq1_read_count = len(fq1)
        fq2_read_count = len(fq2)
    
        fq1_read_total_base = fq1.size
        fq2_read_total_base = fq2.size       
       
        if fq1_read_count!= fq2_read_count:
            print (f'Your input read pair of {sample} have different read counts, you will need to do reads QC before running this software')
            sys.exit()    
        else:
            read_count = fq1_read_count + fq2_read_count
            read_base = fq1_read_total_base + fq2_read_total_base
            sample2read_info[sample] = [read_count, read_base]
            
    elif len(metaG_reads_list) / 2 >= 2 and len(metaG_reads_list) % 2 == 0:
        for i in range(0, len(metaG_reads_list), 2):
            j = i + 1
            sample = Path(metaG_reads_list[i]).stem.rsplit('_', 1)[0]
            
            fq1 = pyfastx.Fastq(metaG_reads_list[i])
            fq2 = pyfastx.Fastq(metaG_reads_list[j])

            fq1_read_count = len(fq1)
            fq2_read_count = len(fq2)
        
            fq1_read_total_base = fq1.size
            fq2_read_total_base = fq2.size       
           
            if fq1_read_count!= fq2_read_count:
                print (f'Your input read pair of {sample} have different read counts, you will need to do reads QC before running this software')
                sys.exit()    
            else:
                read_count = fq1_read_count + fq2_read_count
                read_base = fq1_read_total_base + fq2_read_total_base
                sample2read_info[sample] = [read_count, read_base]  
                
    return sample2read_info           
    
def get_virus_normalized_abundance(mapping_result_dir, virus_raw_abundance, virus_normalized_abundance, sample2read_info, sample2read_info_file):
    # Step 1 Get virus_raw_abundance dict
    virus_raw_abundance_df = pd.read_csv(virus_raw_abundance, sep = '\t', index_col = 0)
    virus_raw_abundance_dict = virus_raw_abundance_df.to_dict() # col => row => value    
    
    virus_normalized_abundance_dict = defaultdict(dict) # sample to gn to cov
    
    for sample in list(virus_raw_abundance_dict.keys()):
        for gn in list(virus_raw_abundance_dict[sample].keys()):
            cov = virus_raw_abundance_dict[sample][gn]
            read_count = sample2read_info[sample][0]
            cov_normalized = cov / (read_count / 100000000) # Normalized by 100M metagenomic reads
            virus_normalized_abundance_dict[sample][gn] = cov_normalized
            
    virus_normalized_abundance_df = pd.DataFrame(virus_normalized_abundance_dict)
    
    virus_normalized_abundance_df['MeanCov'] = virus_normalized_abundance_df.mean(numeric_only=True, axis=1)
    
    for col in virus_normalized_abundance_df.columns:
        col_percent = f'{col}.Percent'
        virus_normalized_abundance_df[col_percent] = (virus_normalized_abundance_df[col] / virus_normalized_abundance_df[col].sum()) * 100
        
    virus_normalized_abundance_df.to_csv(virus_normalized_abundance, sep='\t')

    f = open(sample2read_info_file, 'w')
    f.write('Sample\tRead counts\tRead bases\n')
    for sample in sample2read_info:
        f.write(f'{sample}\t{sample2read_info[sample][0]}\t{sample2read_info[sample][1]}\n')
    f.close()    
        
def parse_vibrant_lytic_and_lysogenic_info(vibrant_outdir, metagenomic_scaffold_stem_name, viral_gn_dir):
    # Step 1 Get scf 2 lytic or lysogenic dict
    lysogenic_fasta_addr = f'{vibrant_outdir}/VIBRANT_phages_{metagenomic_scaffold_stem_name}/{metagenomic_scaffold_stem_name}.phages_lysogenic.fna'
    lytic_fasta_addr = f'{vibrant_outdir}/VIBRANT_phages_{metagenomic_scaffold_stem_name}/{metagenomic_scaffold_stem_name}.phages_lytic.fna'
    
    scf2lytic_or_lyso = {} # scf => 'lytic' or 'lysogenic'
    lysogenic_fasta_seq = store_seq(lysogenic_fasta_addr)
    lysogenic_scf2lyso = {x.replace('>', '', 1):'lysogenic' for x in lysogenic_fasta_seq}
    lytic_fasta_seq = store_seq(lytic_fasta_addr)
    lytic_scf2lytic = {x.replace('>', '', 1):'lytic' for x in lytic_fasta_seq}
    scf2lytic_or_lyso.update(lysogenic_scf2lyso)
    scf2lytic_or_lyso.update(lytic_scf2lytic)
    
    # Step 2 Get gn2long_scfs dict
    gn2long_scfs = {} # gn => [long_scfs]
    all_gn_addrs = glob(f'{viral_gn_dir}/*.fasta')
    for gn_addr in all_gn_addrs:
        gn = Path(gn_addr).stem
        gn_seq = store_seq(gn_addr)
        long_scfs = [x.replace('>', '', 1) for x in gn_seq]
        gn2long_scfs[gn] = long_scfs
        
    # Step 3 Get gn 2 lyso and lytic result
    gn2lyso_lytic_result = {} # gn => 'lytic' or 'lysogenic'
    for gn in gn2long_scfs:
        result = 'lytic' # Default is 'lytic'
        long_scfs = gn2long_scfs[gn]
        for long_scf in long_scfs:
            scf = long_scf.split('__', 1)[1]
            if scf2lytic_or_lyso[scf] == 'lysogenic':
                result = 'lysogenic'
            elif scf2lytic_or_lyso[scf] == 'lytic': 
                continue
            elif scf not in scf2lytic_or_lyso:
                sys.exit('Scaffold {scf} was not parsed correctly')
        gn2lyso_lytic_result[gn] = result        
                
    return gn2lyso_lytic_result            
            
def get_checkv_useful_info(CheckV_quality_summary):
    checkv_table = pd.read_csv(CheckV_quality_summary, sep = '\t', index_col = 0)
    checkv_table = checkv_table[['checkv_quality', 'miuvig_quality', 'completeness', 'completeness_method']]
    checkv_dict = checkv_table.to_dict() # col => row => value
    return checkv_dict
    
def get_viral_gn_size_and_scf_no_and_pro_count(viral_gn_dir):
    gn2size_and_scf_no_and_pro_count = {} # gn => [size, scf_no, pro_count]
    all_gn_addrs = glob(f'{viral_gn_dir}/*.fasta')
    for gn_addr in all_gn_addrs:
        gn = Path(gn_addr).stem
        gn_seq = store_seq(gn_addr)
    
        size = 0
        for header in gn_seq:
            size += len(gn_seq[header])
            
        scf_no = len(gn_seq)
        
        gn_faa_addr = gn_addr.replace('.fasta', '.faa', 1)
        gn_faa_seq = store_seq(gn_faa_addr)
        pro_count = len(gn_faa_seq)
        
        gn2size_and_scf_no_and_pro_count[gn] = [size, scf_no, pro_count]
    return gn2size_and_scf_no_and_pro_count    
    
def get_amg_info_for_vb(vibrant_outdir, metagenomic_scaffold_stem_name, viral_gn_dir):
    gn2long_scf2kos = defaultdict(dict) # gn => long_scf => [kos]
    
    # Step 1 Get gn2long_scfs dict
    gn2long_scfs = {} # gn => [long_scfs]
    all_gn_addrs = glob(f'{viral_gn_dir}/*.fasta')
    for gn_addr in all_gn_addrs:
        gn = Path(gn_addr).stem
        gn_seq = store_seq(gn_addr)
        long_scfs = [x.replace('>', '', 1) for x in gn_seq]
        gn2long_scfs[gn] = long_scfs 
        
    # Step 2 Get scf2kos dict
    scf2kos = defaultdict(list) # scf => [kos]
    with open(f'{vibrant_outdir}/VIBRANT_results_{metagenomic_scaffold_stem_name}/VIBRANT_AMG_individuals_{metagenomic_scaffold_stem_name}.tsv','r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('protein\t'):
                tmp = line.split('\t')
                scf, ko = tmp[1], tmp[2]
                scf2kos[scf].append(ko)
                
    # Step 3 Get gn2long_scf2kos dict
    for gn in gn2long_scfs:
        for long_scf in gn2long_scfs[gn]:
            scf = long_scf.split('__', 1)[1]
            kos = scf2kos[scf]
            gn2long_scf2kos[gn][long_scf] = kos
            
    return gn2long_scf2kos  

def get_amg_info_for_vs_and_dvf(args, viral_gn_dir):
    gn2long_scf2kos = defaultdict(dict) # gn => long_scf => [kos]
    
    # Step 1 Get gn2long_scfs dict
    gn2long_scfs = {} # gn => [long_scfs]
    all_gn_addrs = glob(f'{viral_gn_dir}/*.fasta')
    for gn_addr in all_gn_addrs:
        gn = Path(gn_addr).stem
        gn_seq = store_seq(gn_addr)
        long_scfs = [x.replace('>', '', 1) for x in gn_seq]
        gn2long_scfs[gn] = long_scfs 
        
    # Step 2 Get scf2kos dict
    scf2kos = defaultdict(list) # scf => [kos]
    annotation_file = ''
    if args['identify_method'] == 'vs':
        annotation_file = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.annotation.txt')
    elif args['identify_method'] == 'dvf':
        annotation_file = os.path.join(args['dvf_outdir'], 'final_dvf_virus.annotation.txt')
    with open(annotation_file ,'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('protein\t'):
                tmp = line.split('\t')
                scf, ko = tmp[1], tmp[2]
                scf2kos[scf].append(ko)
                
    # Step 3 Get gn2long_scf2kos dict
    for gn in gn2long_scfs:
        for long_scf in gn2long_scfs[gn]:
            scf = long_scf.split('__', 1)[1]
            kos = scf2kos[scf]
            gn2long_scf2kos[gn][long_scf] = kos
            
    return gn2long_scf2kos     
            
def get_amg_statics(gn2long_scf2kos):
    gn2amg_statics = {} # gn => amg_statics; for example, K00018(3);K01953(4)
    for gn in gn2long_scf2kos:
        amg_statics_list = []
        ko2hit_num = {} # ko => hit_num
        for long_scf in gn2long_scf2kos[gn]:
            kos = gn2long_scf2kos[gn][long_scf]
            for ko in kos:
                ko2hit_num[ko] = ko2hit_num.get(ko, 0) + 1
        for ko in ko2hit_num:
            hit_num = ko2hit_num[ko]
            amg_statics_list.append(f'{ko}({hit_num})')
        gn2amg_statics[gn] = ';'.join(amg_statics_list)  
        
    return gn2amg_statics    
    
def get_virus_summary_info(checkv_dict, gn2lyso_lytic_result, gn2size_and_scf_no_and_pro_count, gn2amg_statics, virus_summary_info):
    gns = list(gn2size_and_scf_no_and_pro_count.keys())
    
    virus_summary_info_dict = defaultdict(dict) # parameter => gn => value
    virus_summary_info_dict.update(checkv_dict)
    
    for gn in gns:
        lytic_state = gn2lyso_lytic_result.get(gn, '')
        virus_summary_info_dict['lytic_state'][gn] = lytic_state
        virus_summary_info_dict['genome_size'][gn] = gn2size_and_scf_no_and_pro_count[gn][0]
        virus_summary_info_dict['scaffold_num'][gn] = gn2size_and_scf_no_and_pro_count[gn][1]
        virus_summary_info_dict['protein_count'][gn] = gn2size_and_scf_no_and_pro_count[gn][2]
        virus_summary_info_dict['AMG_KOs'][gn] = gn2amg_statics[gn]
        
    virus_summary_info_df = pd.DataFrame(virus_summary_info_dict) 

    virus_summary_info_df = virus_summary_info_df[['genome_size', 'scaffold_num', 'protein_count', 'AMG_KOs', 'lytic_state', 'checkv_quality', 'miuvig_quality', 'completeness', 'completeness_method']]
    virus_summary_info_df.to_csv(virus_summary_info, sep='\t')

def get_run_input_arguments(args):
    command = f"{os.path.join(args['root_dir'], 'ViWrap')} run "
    argu_items = []
    if args['input_metagenome'] != 'none': argu_items.append('--input_metagenome' + ' ' + args['input_metagenome'])
    if args['input_reads'] != 'none': argu_items.append('--input_reads' + ' ' + args['input_reads'])  
    argu_items.append('--out_dir' + ' ' + args['out_dir'])
    argu_items.append('--db_dir' + ' ' + args['db_dir'])
    argu_items.append('--identify_method' + ' ' + args['identify_method'])
    if args['conda_env_dir'] != 'none': argu_items.append('--conda_env_dir' + ' ' + args['conda_env_dir'])
    argu_items.append('--threads' + ' ' + args['threads'])
    if args['virome']: argu_items.append('--virome')
    argu_items.append('--input_length_limit' + ' ' + args['input_length_limit'])
    if args['custom_MAGs_dir'] != 'none': argu_items.append('--custom_MAGs_dir' + ' ' + args['custom_MAGs_dir'])
    
    command += " ".join(argu_items)
    return command
    
def combine_iphop_results(args, combined_host_pred_to_genome_result, combined_host_pred_to_genus_result):
    host_pred_to_genome_m90 = os.path.join(args['iphop_outdir'], "Host_prediction_to_genome_m90.csv")
    host_pred_to_genus_m90 = os.path.join(args['iphop_outdir'], "Host_prediction_to_genus_m90.csv")
    
    host_pred_to_genome_header = ''
    host_pred_to_genome_m90_result = [] # Store each line
    with open(host_pred_to_genome_m90, 'r') as lines:
        for line in lines:
            line = line.strip('\n')
            if line.startswith('Virus,'):
                host_pred_to_genome_header = line
            else:
                host_pred_to_genome_m90_result.append(line)
    lines.close()            

    host_pred_to_genus_header = ''
    host_pred_to_genus_m90_result = [] # Store each line
    with open(host_pred_to_genus_m90, 'r') as lines:
        for line in lines:
            line = line.strip('\n')
            if line.startswith('Virus,'):
                host_pred_to_genus_header = line
            else:
                host_pred_to_genus_m90_result.append(line)
    lines.close()                

    if args['custom_MAGs_dir'] != 'none':
        host_pred_to_genome_m90_custom = os.path.join(args['iphop_custom_outdir'], "Host_prediction_to_genome_m90.csv")
        host_pred_to_genus_m90_custom = os.path.join(args['iphop_custom_outdir'], "Host_prediction_to_genus_m90.csv")
        
        with open(host_pred_to_genome_m90_custom, 'r') as lines:
            for line in lines:
                line = line.strip('\n')
                if not line.startswith('Virus,'):
                    if line not in host_pred_to_genome_m90_result:
                        host_pred_to_genome_m90_result.append(line)
        lines.close()            
                    
        with open(host_pred_to_genus_m90_custom, 'r') as lines:
            for line in lines:
                line = line.strip('\n')
                if not line.startswith('Virus,'):
                    if line not in host_pred_to_genus_m90_result:
                        host_pred_to_genus_m90_result.append(line)
        lines.close()  

    f = open(combined_host_pred_to_genome_result, 'w')
    f.write(host_pred_to_genome_header + '\n')
    for line in host_pred_to_genome_m90_result:
        f.write(line + '\n')
    f.close()

    f = open(combined_host_pred_to_genus_result , 'w')
    f.write(host_pred_to_genus_header + '\n')
    for line in host_pred_to_genus_m90_result:
        f.write(line + '\n')
    f.close()  

def get_virus_genome_annotation_result(args):
    if args['identify_method'] == 'vb':
        # Step 1 Store annotation result
        vibrant_annotation_result_file = os.path.join(args['vibrant_outdir'],f"VIBRANT_results_{Path(args['input_metagenome']).stem}",f"VIBRANT_annotations_{Path(args['input_metagenome']).stem}.tsv")
        
        vibrant_annotation_result = {} # protein => [items in each line]
        vibrant_annotation_result_header = ''
        with open (vibrant_annotation_result_file, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                if line.startswith('protein\t'):
                    vibrant_annotation_result_header = line
                else:
                    tmp = line.split('\t')
                    protein = tmp[0]
                    vibrant_annotation_result[protein] = tmp
        lines.close()            
                    
        # Step 2 Store gn2long_proteins and long_protein2gn dict
        gn2long_proteins = {} # gn => [long_proteins]
        long_protein2gn = {} # long_protein => gn
        
        all_faa_addrs = glob(os.path.join(args['viwrap_summary_outdir'],'Virus_genomes_files/*.faa'))
        for faa_addr in all_faa_addrs:
            gn = Path(faa_addr).stem
            faa_seqs = store_seq(faa_addr)
            long_proteins = [x.replace('>', '', 1).split('\t', 1)[0] for x in faa_seqs]
            gn2long_proteins[gn] = long_proteins
            
            for long_protein in long_proteins:
                long_protein2gn[long_protein] = gn
                
        # Step 3 Get new VIBRANT annotation result
        vibrant_annotation_result_new = {} # long_protein => [items in each line]
        for long_protein in long_protein2gn:
            protein = long_protein.split('__', 1)[1]
            items = vibrant_annotation_result[protein]
            items[0] = long_protein
            items[1] = long_protein.rsplit('_', 1)[0]
            items.insert(0, long_protein2gn[long_protein])
            vibrant_annotation_result_new[long_protein] = items
            
        # Step 4 Write down the result
        result = os.path.join(args['viwrap_summary_outdir'],'Virus_annotation_results.txt')
        f = open(result, 'w')
        f.write('viral genome\t' + vibrant_annotation_result_header + '\n')
        for long_protein in vibrant_annotation_result_new:
            line = '\t'.join(vibrant_annotation_result_new[long_protein])
            f.write(line + '\n')
        f.close() 
    elif args['identify_method'] == 'vs' or args['identify_method'] == 'dvf':
        # Step 1 Store annotation result
        annotation_result_file = ''
        if args['identify_method'] == 'vs':
            annotation_result_file = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.annotation.txt')
        elif args['identify_method'] == 'dvf': 
            annotation_result_file = os.path.join(args['dvf_outdir'], 'final_dvf_virus.annotation.txt')
        
        annotation_result = {} # protein => [items in each line]
        annotation_result_header = ''
        with open (annotation_result_file, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                if line.startswith('protein\t'):
                    annotation_result_header = line
                else:
                    tmp = line.split('\t')
                    protein = tmp[0]
                    annotation_result[protein] = tmp
        lines.close() 

        # Step 2 Store gn2long_proteins and long_protein2gn dict
        gn2long_proteins = {} # gn => [long_proteins]
        long_protein2gn = {} # long_protein => gn
        
        all_faa_addrs = glob(os.path.join(args['viwrap_summary_outdir'],'Virus_genomes_files/*.faa'))
        for faa_addr in all_faa_addrs:
            gn = Path(faa_addr).stem
            faa_seqs = store_seq(faa_addr)
            long_proteins = [x.replace('>', '', 1) for x in faa_seqs]
            gn2long_proteins[gn] = long_proteins
            
            for long_protein in long_proteins:
                long_protein2gn[long_protein] = gn 

        # Step 3 Get new VIBRANT annotation result
        annotation_result_new = {} # long_protein => [items in each line]
        for long_protein in long_protein2gn:
            protein = long_protein.split('__', 1)[1]
            items = annotation_result[protein]
            items[0] = long_protein
            items[1] = long_protein.rsplit('_', 1)[0]
            items.insert(0, long_protein2gn[long_protein])
            annotation_result_new[long_protein] = items
            
        # Step 4 Write down the result
        result = os.path.join(args['viwrap_summary_outdir'],'Virus_annotation_results.txt')
        f = open(result, 'w')
        f.write('viral genome\t' + annotation_result_header + '\n')
        for long_protein in annotation_result_new:
            line = '\t'.join(annotation_result_new[long_protein])
            f.write(line + '\n')
        f.close()                      

def screen_virsorter2_result(args, keep1_list_file, keep2_list_file, discard_list_file, manual_check_list_file):
    seq2info = defaultdict(list) # seq => [length, score, hallmark, viral_gene, host_gene]
    # Step 1 Parse final_viral_score file from VirSorter2 pass2 folder
    final_viral_score = os.path.join(args['virsorter_outdir'], 'pass2/final-viral-score.tsv')
    with open(final_viral_score, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('seqname\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark = tmp[0], float(tmp[5]), float(tmp[3]), float(tmp[6])
                seq2info[seq].append(length)
                seq2info[seq].append(score)
                seq2info[seq].append(hallmark)
    lines.close()

    # Step 2 Parse quality_summary file from VirSorter2 CheckV_result_2nd folder
    quality_summary = os.path.join(args['virsorter_outdir'], 'CheckV_result_2nd/quality_summary.tsv')
    with open(quality_summary, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('contig_id\t'):
                tmp = line.split('\t')
                seq, viral_gene, host_gene = tmp[0], float(tmp[5]), float(tmp[6])
                seq2info[seq].append(viral_gene)
                seq2info[seq].append(host_gene)
    lines.close() 

    # Step 3 Make keep_list_file, discard_list_file, manual_check_list_file
    keep1_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    keep2_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    discard_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    manual_check_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    
    for seq in seq2info:
        length, score, hallmark, viral_gene, host_gene = seq2info[seq][0], seq2info[seq][1], seq2info[seq][2], seq2info[seq][3], seq2info[seq][4] 
        if viral_gene > 0:
            keep1_list[seq] = [length, score, hallmark, viral_gene, host_gene]
        elif viral_gene == 0 and (score >= 0.95 or hallmark > 2 or host_gene == 0):
            keep2_list[seq] = [length, score, hallmark, viral_gene, host_gene]
        elif viral_gene == 0 and host_gene == 1 and length >= 10000:
            manual_check_list[seq] = [length, score, hallmark, viral_gene, host_gene]
            
    for seq in manual_check_list:
        if seq in keep1_list:
            del keep1_list[seq]
        if seq in keep2_list:
            del keep2_list[seq]    
            
    for seq in seq2info:
        if seq not in keep1_list and seq not in keep2_list and seq not in manual_check_list:
            discard_list[seq] = seq2info[seq]
        
            
    f = open(keep1_list_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in keep1_list:
        f.write(seq + '\t' + '\t'.join(str(item) for item in keep1_list[seq]) + '\n')
    f.close()  
    
    f = open(keep2_list_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in keep2_list:
        f.write(seq + '\t' + '\t'.join(str(item) for item in keep2_list[seq]) + '\n')
    f.close()      

    f = open(discard_list_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in discard_list:
        f.write(seq + '\t' + '\t'.join(str(item) for item in discard_list[seq]) + '\n')
    f.close() 

    f = open(manual_check_list_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in manual_check_list:
        f.write(seq + '\t' + '\t'.join(str(item) for item in manual_check_list[seq]) + '\n')
    f.close()     
    
def get_keep2_mc_seq(args, keep2_list_file, manual_check_list_file, keep2_fasta, manual_check_fasta):
    # Step 1 Store keep2_list, manual_check_list
    keep2_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    with open(keep2_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                keep2_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()

    manual_check_list = {} # seq => [length, score, hallmark, viral_gene, host_gene] 
    with open(manual_check_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                manual_check_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()  

    # Step 2 Make keep2_fasta, manual_check_fasta
    all_seq = store_seq(os.path.join(args['virsorter_outdir'], 'pass2/final-viral-combined.fa'))
    
    all_seq_keep2 = {}
    for header in all_seq:
        header_wo_array = header.replace('>', '', 1)
        if header_wo_array in keep2_list:
            all_seq_keep2[header] = all_seq[header]
            
    all_seq_manual_check = {}
    for header in all_seq:
        header_wo_array = header.replace('>', '', 1)
        if header_wo_array in manual_check_list:
            all_seq_manual_check[header] = all_seq[header] 
            
    write_down_seq(all_seq_keep2, keep2_fasta) 
    write_down_seq(all_seq_manual_check, manual_check_fasta) 

def get_keep2_vb_passed_list(args, keep2_vb_result, keep2_list_vb_passed_file):
    # Step 1 Store keep2_list
    keep2_list_file = os.path.join(args['virsorter_outdir'], 'keep2_list.txt')
    
    keep2_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    with open(keep2_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                keep2_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()
  
    keep2_list_vb_passed = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    keep2_vb_result_seq = store_seq(keep2_vb_result)
    for header in keep2_vb_result_seq:
        header_wo_array = header.replace('>', '', 1)
        if '_fragment' in header_wo_array:
            header_wo_array = header_wo_array.rsplit('_fragment', 1)[0]
        keep2_list_vb_passed[header_wo_array] = keep2_list[header_wo_array]
        
    f = open(keep2_list_vb_passed_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in keep2_list_vb_passed:
        f.write(seq + '\t' + '\t'.join(str(item) for item in keep2_list_vb_passed[seq]) + '\n')
    f.close()  

def get_manual_check_vb_passed_list(args, manual_check_vb_result, manual_check_list_vb_passed_file):
    # Step 1 Store manual_check_list
    manual_check_list_file = os.path.join(args['virsorter_outdir'], 'manual_check_list.txt')
    
    manual_check_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    with open(manual_check_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                manual_check_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()
  
    manual_check_list_vb_passed = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    manual_check_vb_result_seq = store_seq(manual_check_vb_result)
    for header in manual_check_vb_result_seq:
        header_wo_array = header.replace('>', '', 1)
        if '_fragment' in header_wo_array:
            header_wo_array = header_wo_array.rsplit('_fragment', 1)[0]
        manual_check_list_vb_passed[header_wo_array] = manual_check_list[header_wo_array]
        
    f = open(manual_check_list_vb_passed_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in manual_check_list_vb_passed:
        f.write(seq + '\t' + '\t'.join(str(item) for item in manual_check_list_vb_passed[seq]) + '\n')
    f.close()     
    
def get_final_vs2_virus(args, keep1_list_file, keep2_list_vb_passed_file, manual_check_list_vb_passed_file, final_vs2_virus_fasta_file):
    # Step 1 Store keep1_list,  keep2_list_vb_passed, manual_check_list_vb_passed
    keep1_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    with open(keep1_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                keep1_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()     
    
    keep2_list_vb_passed = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    if os.path.exists(keep2_list_vb_passed_file):
        with open(keep2_list_vb_passed_file, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                if not line.startswith('#seq\t'):
                    tmp = line.split('\t')
                    seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                    keep2_list_vb_passed[seq] = [length, score, hallmark, viral_gene, host_gene]
        lines.close() 

    manual_check_list_vb_passed = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    if os.path.exists(manual_check_list_vb_passed_file):    
        with open(manual_check_list_vb_passed_file, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                if not line.startswith('#seq\t'):
                    tmp = line.split('\t')
                    seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                    manual_check_list_vb_passed[seq] = [length, score, hallmark, viral_gene, host_gene]
        lines.close()  

    # Step 2 Make final_vs2_virus.fasta
    all_seq = store_seq(os.path.join(args['virsorter_outdir'], 'pass2/final-viral-combined.fa'))
    
    all_seq_final = {}
    for header in all_seq:
        header_wo_array = header.replace('>', '', 1)
        if header_wo_array in keep1_list or header_wo_array in keep2_list_vb_passed or header_wo_array in manual_check_list_vb_passed:
            all_seq_final[header] = all_seq[header] 
            
    write_down_seq(all_seq_final, final_vs2_virus_fasta_file)    
    
def get_dvf_result_seq(args, final_dvf_virus_fasta_file):
    # Step 1 Store and filter dvfpred.txt
    dvf_passed_seq = [] 
    with open(os.path.join(args['dvf_outdir'], f"{Path(args['input_metagenome']).stem}.fasta_gt{args['input_length_limit']}bp_dvfpred.txt"),'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('name\t'):
                tmp = line.split('\t')
                seq = tmp[0]
                score = tmp[2]
                pvalue = tmp[3]
                if float(score) >= 0.95 and float(pvalue) < 0.05: 
                    dvf_passed_seq.append(seq)
                else:
                    continue
             
    # Step 2 get the final_dvf_virus_fasta_file
    all_seq = store_seq(args['input_metagenome'])
    
    all_seq_final = {}
    for header in all_seq:
        header_wo_array = header.replace('>', '', 1)
        if header_wo_array in dvf_passed_seq:
            all_seq_final[header] = all_seq[header] 
            
    write_down_seq(all_seq_final, final_dvf_virus_fasta_file)  
                 
    
    
            
    
    
                
    
    
                
    
    
    
    
    
            
                    
                    
        
        
    
    
    
        
        
        
        
    
            
    
            
    
    
    
    