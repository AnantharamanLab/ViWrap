#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    warnings.filterwarnings("ignore")
    from pathlib import Path
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call
    from Bio import SeqIO
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

def write_down_seq(seq_dict, path_to_file): 
    # Two inputs are required:
    # (1) The dict of the sequence
    # (2) The path that you want to write your sequence down
    
    seq_file = open(path_to_file,"w")
    for head in seq_dict:
        seq_file.write(head + "\n")
        seq_file.write(seq_dict[head] + "\n")
    seq_file.close()
    
def dl_refseq_viral_protein(tax_classification_db_dir):
    dl_cmd = []
    for i in range(1, 4):
        each_dl_cmd = f'wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{i}.protein.faa.gz -O {tax_classification_db_dir}/viral.{i}.protein.faa.gz'
        each_gzip_cmd = f'gzip -d {tax_classification_db_dir}/viral.{i}.protein.faa.gz'
        each_cmd = each_dl_cmd + ";" + each_gzip_cmd
        dl_cmd.append(each_cmd)
    
    n = 4 # The number of parallel processes
    for j in range(max(int(len(dl_cmd)/n + 1), 1)):
        procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in dl_cmd[j*n: min((j+1)*n, len(dl_cmd))] ]
        for p in procs:
            p.wait()  
            
    combind_cmd = f'cat {tax_classification_db_dir}/viral.*.protein.faa > {tax_classification_db_dir}/NCBI_RefSeq_viral.faa'
    rm_cmd = f'rm {tax_classification_db_dir}/viral.*.protein.faa'
    os.system(combind_cmd + ";" + rm_cmd)
    
def dl_refseq_viral_protein_gpff(tax_classification_db_dir):
    dl_cmd = []
    #os.mkdir(tax_classification_db_dir)
    for i in range(1, 4):
        each_dl_cmd = f'wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{i}.protein.gpff.gz -O {tax_classification_db_dir}/viral.{i}.protein.gpff.gz'
        each_gzip_cmd = f'gzip -d {tax_classification_db_dir}/viral.{i}.protein.gpff.gz'
        each_cmd = each_dl_cmd + ";" + each_gzip_cmd
        dl_cmd.append(each_cmd)
    
    n = 4 # The number of parallel processes
    for j in range(max(int(len(dl_cmd)/n + 1), 1)):
        procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in dl_cmd[j*n: min((j+1)*n, len(dl_cmd))] ]
        for p in procs:
            p.wait()  
            
    combind_cmd = f'cat {tax_classification_db_dir}/viral.*.protein.gpff > {tax_classification_db_dir}/NCBI_RefSeq_viral.gpff'
    rm_cmd = f'rm {tax_classification_db_dir}/viral.*.protein.gpff'
    os.system(combind_cmd + ";" + rm_cmd)
    
def parse_gpff(tax_classification_db_dir):
    pro2tax = {}
    gpff = f'{tax_classification_db_dir}/NCBI_RefSeq_viral.gpff'
    f = open(gpff, "r")
    for gp_record in SeqIO.parse(f, 'genbank'):
        acc = gp_record.annotations['accessions'][0]
        organism = gp_record.annotations['organism']
        tax_line = (";").join(gp_record.annotations['taxonomy'])
        pro2tax[acc] = tax_line + ";" + organism
        
    f.close()

    fo = open(f'{tax_classification_db_dir}/NCBI_RefSeq_viral_protein2NCBI_tax.txt', "w")
    for pro in pro2tax:
        tax = pro2tax[pro]
        fo.write(f'{pro}\t{tax}\n')
    fo.close()    
    
def grep_NCBI_RefSeq_viral_proteins_w_tax(tax_classification_db_dir):
    pro_seq = store_seq(f'{tax_classification_db_dir}/NCBI_RefSeq_viral.faa')
    accessions_w_tax = set()
    
    with open(f'{tax_classification_db_dir}/NCBI_RefSeq_viral_protein2NCBI_tax.txt',"r") as lines:
        for line in lines:
            line = line.rstrip("\n")
            acc = line.split("\t", 1)[0]
            accessions_w_tax.add(acc)
    lines.close()        
    
    pro_seq_w_tax = {}
    for pro in pro_seq:
        acc = pro.replace(">", "").split(".", 1)[0]
        if acc in accessions_w_tax:
            acc_w_arrow = ">" + acc
            pro_seq_w_tax[acc_w_arrow] = pro_seq[pro]
    
    write_down_seq(pro_seq_w_tax, f'{tax_classification_db_dir}/NCBI_RefSeq_viral.faa')

def reformat_NCBI_tax_to_ICTV_8_rank_tax(tax_classification_db_dir, ictv_tax_info, pro2ictv_8_rank_tax):
    # Step 1 Store NCBI tax and dereplicate it
    ID2NCBI_tax = {} # ID => NCBI_tax
    NCBI_tax_dict = {} # NCBI_tax => species
    
    with open(f'{tax_classification_db_dir}/NCBI_RefSeq_viral_protein2NCBI_tax.txt', "r") as lines:
        for line in lines:
            line = line.rstrip("\n")
            ID = line.split("\t")[0]
            NCBI_tax = line.split("\t")[1]
            species = NCBI_tax.rsplit(";", 1)[1]
            ID2NCBI_tax[ID] = NCBI_tax
            NCBI_tax_dict[NCBI_tax] = species
            
    # Step 2 Store ICTV tax
    Realm_set = Subrealm_set = Kingdom_set = Subkingdom_set = set()
    Phylum_set = Subphylum_set = Class_set = Subclass_set = set()
    Order_set = Suborder_set = Family_set = Subfamily_set =  set()
    Genus_set = Subgenus_set = Species_set = set()

    String2rank = {} # string => rank (i.e., realm, subrealm)
    
    with open(ictv_tax_info, "r") as lines:
        for line in lines:
            line = line.rstrip("\n")
            if not line.startswith("Sort") and not line.startswith("#"):
                tmp = line.split("\t")
                sort_id, realm, subrealm, kingdom, subkingdom = tmp[0], tmp[1], tmp[2], tmp[3], tmp[4]
                phylum, subphylum, class_, subclass, order, suborder = tmp[5], tmp[6], tmp[7], tmp[8], tmp[9], tmp[10]
                family, subfamily, genus, subgenus, species, genome_composition = tmp[11], tmp[12], tmp[13], tmp[14], tmp[15], tmp[16]
                
                Realm_set.add(realm)
                Subrealm_set.add(subrealm)
                Kingdom_set.add(kingdom)
                Subkingdom_set.add(subkingdom)
                Phylum_set.add(phylum)
                Subphylum_set.add(subphylum)
                Class_set.add(class_)
                Subclass_set.add(subclass)
                Order_set.add(order)
                Suborder_set.add(suborder)
                Family_set.add(family)
                Subfamily_set.add(subfamily)
                Genus_set.add(genus)
                Subgenus_set.add(subgenus)
                Species_set.add(species)
                
                String2rank[realm] = 'Realm'
                String2rank[subrealm] = 'Subrealm'
                String2rank[kingdom] = 'Kingdom'
                String2rank[subkingdom] = 'Subkingdom'
                String2rank[phylum] = 'Phylum'
                String2rank[subphylum] = 'Subphylum'
                String2rank[class_] = 'Class'
                String2rank[subclass] = 'Subclass'
                String2rank[order] = 'Order'
                String2rank[suborder] = 'Suborder'
                String2rank[family] = 'Family'
                String2rank[subfamily] = 'Subfamily'
                String2rank[genus] = 'Genus'
                String2rank[subgenus] = 'Subgenus'
                String2rank[species] = 'Species'
    lines.close()
                
    # Step 3. Reformat the NCBI tax to 8-rank ICTV tax
    NCBI_tax2ictv_8_rank_tax = {} # NCBI_tax => ictv_8_rank_tax
    
    for NCBI_tax in NCBI_tax_dict:
        species = NCBI_tax_dict[NCBI_tax] # The "species" is according to the NCBI species
        NCBI_tax = NCBI_tax.replace("Viruses;", "", 1) # Delete the "Viruses;" in the front
        ictv_8_rank_tax = ''
        realm = kingdom = phylum = class_ = order = family = genus = ''
        
        NCBI_tax_list = NCBI_tax.split(";")
        for i in range((len(NCBI_tax_list) - 1)):
            if NCBI_tax_list[i] in String2rank:
                if String2rank[NCBI_tax_list[i]] == 'Realm':
                    realm = NCBI_tax_list[i]
                elif String2rank[NCBI_tax_list[i]] == 'Kingdom':
                    kingdom = NCBI_tax_list[i]
                elif String2rank[NCBI_tax_list[i]] == 'Phylum':
                    phylum = NCBI_tax_list[i]
                elif String2rank[NCBI_tax_list[i]] == 'Class':     
                    class_ = NCBI_tax_list[i]
                elif String2rank[NCBI_tax_list[i]] == 'Order':
                    order = NCBI_tax_list[i]
                elif String2rank[NCBI_tax_list[i]] == 'Family':
                    family = NCBI_tax_list[i]
                elif String2rank[NCBI_tax_list[i]] == 'Genus':     
                    genus = NCBI_tax_list[i]                    
        ictv_8_rank_tax = f'{realm};{kingdom};{phylum};{class_};{order};{family};{genus};{species}'
        NCBI_tax2ictv_8_rank_tax[NCBI_tax] = ictv_8_rank_tax # Note that here NCBI_tax does not contain "Viruses;" in the front
        
    # Step 4. Write down result
    f = open(pro2ictv_8_rank_tax, "w")
    for ID in ID2NCBI_tax:
        NCBI_tax = ID2NCBI_tax[ID]
        NCBI_tax = NCBI_tax.replace("Viruses;", "", 1)
        ictv_8_rank_tax = NCBI_tax2ictv_8_rank_tax[NCBI_tax]
        f.write(f'{ID}\t{ictv_8_rank_tax}\n')
    f.close()   

def make_diamond_db(tax_classification_db_dir):
    faa = f'{tax_classification_db_dir}/NCBI_RefSeq_viral.faa'
    faa_stem = faa.rsplit(".", 1)[0]
    cmd = f'diamond makedb --in {faa} --db {faa_stem}'
    os.system(cmd)

def remove(tax_classification_db_dir):    
    remove_cmds = []
    remove_cmds.append(f'rm {tax_classification_db_dir}/NCBI_RefSeq_viral.gpff')
    remove_cmds.append(f'rm {tax_classification_db_dir}/NCBI_RefSeq_viral_protein2NCBI_tax.txt')
    os.system(';'.join(remove_cmds))
    
def get_vog_marker_table(vog_marker_table):
    vog_marker_list = {} # vog => tax
    with open(vog_marker_table, "r") as lines:
        for line in lines:
            line = line.rstrip("\n")
            if not line.startswith("#"):
                vog, tax = line.split("\t")[0], line.split("\t")[2]
                vog_marker_list[vog] = tax
    lines.close()
    return vog_marker_list
    
def get_marker_vog_hmm(vog_marker_list, tax_classification_db_dir):
    # Step 1 Download whole VOG HMMs (VOG 97)
    wget_cmd = f'wget http://fileshare.csb.univie.ac.at/vog/vog97/vog.hmm.tar.gz -O {tax_classification_db_dir}/vog.hmm.tar.gz'
    os.system(wget_cmd)
    # Step 2 Pick marker VOG HMMs, concatenate, and press
    os.mkdir(f'{tax_classification_db_dir}/tmp')
    unzip_cmd = f'tar xzf {tax_classification_db_dir}/vog.hmm.tar.gz --directory {tax_classification_db_dir}/tmp'
    os.system(unzip_cmd)
    
    marker_hmms = []
    for vog in vog_marker_list:
        marker_hmms.append(f'{tax_classification_db_dir}/tmp/{vog}.hmm')
    
    cat_cmd = f'cat {" ".join(marker_hmms)} > {tax_classification_db_dir}/marker_VOG.hmm'
    os.system(cat_cmd)
    
    press_cmd = f'hmmpress {tax_classification_db_dir}/marker_VOG.hmm'
    os.system(press_cmd)
    
    os.system(f'rm -rf {tax_classification_db_dir}/tmp')
    os.system(f'rm {tax_classification_db_dir}/vog.hmm.tar.gz')
    