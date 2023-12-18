#!/usr/bin/env python3

try:
    import warnings
    import sys
    from pathlib import Path
    import os
    import re
    warnings.filterwarnings("ignore")
    from glob import glob
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 


def store_seq(input_seq_file):
    head = "" 
    seq_dict = {}
    
    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.rstrip("\n")
            if ">" in line:
                if " " in line or "\t" in line:
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
            
    return seq_dict  

def get_LCA_from_taxs(taxs):
    Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species = {}, {}, {}, {}, {}, {}, {}, {}
    for tax in taxs:
        tmp = tax.split(';')
        Realm[tmp[0]] = Realm.get(tmp[0], 0) + 1
        Kingdom[tmp[1]] = Kingdom.get(tmp[1], 0) + 1
        Phylum[tmp[2]] = Phylum.get(tmp[2], 0) + 1
        Class[tmp[3]] = Class.get(tmp[3], 0) + 1
        Order[tmp[4]] = Order.get(tmp[4], 0) + 1
        Family[tmp[5]] = Family.get(tmp[5], 0) + 1
        Genus[tmp[6]] = Genus.get(tmp[6], 0) + 1
        Species[tmp[7]] = Species.get(tmp[7], 0) + 1
    
    lca_list = []
    if len(Realm) == 1 and len(Kingdom) != 1:
        lca_list = [list(Realm.keys())[0]]
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) != 1:
        lca_list = [list(Realm.keys())[0], list(Kingdom.keys())[0]]
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) != 1:
        lca_list = [list(Realm.keys())[0], list(Kingdom.keys())[0], list(Phylum.keys())[0]]
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) != 1:
        lca_list = [list(Realm.keys())[0], list(Kingdom.keys())[0], list(Phylum.keys())[0], list(Class.keys())[0]]
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) == 1 and len(Family) != 1:
        lca_list = [list(Realm.keys())[0], list(Kingdom.keys())[0], list(Phylum.keys())[0], list(Class.keys())[0], list(Order.keys())[0]]
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) == 1 and len(Family) == 1 and len(Genus) != 1:
        lca_list = [list(Realm.keys())[0], list(Kingdom.keys())[0], list(Phylum.keys())[0], list(Class.keys())[0], list(Order.keys())[0], list(Family.keys())[0]]
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) == 1 and len(Family) == 1 and len(Genus) == 1 and len(Species) != 1:
        lca_list = [list(Realm.keys())[0], list(Kingdom.keys())[0], list(Phylum.keys())[0], list(Class.keys())[0], list(Order.keys())[0], list(Family.keys())[0], list(Genus.keys())[0]]
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) == 1 and len(Family) == 1 and len(Genus) == 1 and len(Species) == 1:
        lca_list = [list(Realm.keys())[0], list(Kingdom.keys())[0], list(Phylum.keys())[0], list(Class.keys())[0], list(Order.keys())[0], list(Family.keys())[0], list(Genus.keys())[0], list(Species.keys())[0]]

    lca_list_full = 'NA;NA;NA;NA;NA;NA;NA;NA'.split(';')
    if lca_list:
        for i in range(len(lca_list)):
            lca_list_full[i] = lca_list[i]
            
    return ';'.join(lca_list_full)     
     
def get_tax_from_genomad_result(genomad_summary_file, vRhyme_best_bin_dir_modified, vRhyme_unbinned_viral_gn_dir, tax_genomad_tax_output):
    viral_scf2tax = {}
    with open(genomad_summary_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('seq_name'):
                tmp = line.split('\t')
                viral_scf, tax = tmp[0], tmp[-1]
                if tax != 'Unclassified':
                    tax_list = tax.split(';')
                    tax_list.pop(0)
                    while len(tax_list) < 8:
                        tax_list.append("NA")
                    tax_new = ';'.join(tax_list)
                    viral_scf2tax[viral_scf] = tax_new
                    
    gn2viral_scf_list = {}
    vRhyme_best_bin_addrs = glob(f"{vRhyme_best_bin_dir_modified}/*.fasta")
    for each_vRhyme_best_bin_addr in vRhyme_best_bin_addrs:
        each_vRhyme_best_bin_seq = store_seq(each_vRhyme_best_bin_addr)
        gn = Path(each_vRhyme_best_bin_addr).stem
        if 'vRhyme' in gn:  # If the gn is a vRhyme bin
            gn2viral_scf_list[gn] = [x.replace('>', '', 1).split('__', 1)[1] for x in each_vRhyme_best_bin_seq]
        else: # If the gn is not a vRhyme bin, but a single-scaffold virus
            gn2viral_scf_list[gn] = [x.replace('>', '', 1) for x in each_vRhyme_best_bin_seq]    
    
    vRhyme_unbinned_viral_gn_addrs = glob(f"{vRhyme_unbinned_viral_gn_dir}/*.fasta")
    for each_vRhyme_unbinned_viral_gn_addr in vRhyme_unbinned_viral_gn_addrs:
        each_vRhyme_unbinned_viral_gn_seq = store_seq(each_vRhyme_unbinned_viral_gn_addr)
        gn = Path(each_vRhyme_unbinned_viral_gn_addr).stem
        if 'vRhyme' in gn:  # If the gn is a vRhyme bin  
            gn2viral_scf_list[gn] = [x.replace('>', '', 1).split('__', 1)[1] for x in each_vRhyme_unbinned_viral_gn_seq]
        else: # If the gn is not a vRhyme bin, but a single-scaffold virus
            gn2viral_scf_list[gn] = [x.replace('>', '', 1) for x in each_vRhyme_unbinned_viral_gn_seq]
                                
    gn2lca_tax = {}
    for gn in gn2viral_scf_list:
        taxs = []
        for viral_scf in gn2viral_scf_list[gn]:
            if viral_scf in viral_scf2tax:
                tax = viral_scf2tax[viral_scf]
                taxs.append(tax)
        lca_tax = get_LCA_from_taxs(taxs)           
        if lca_tax != 'NA;NA;NA;NA;NA;NA;NA;NA':
            gn2lca_tax[gn] = lca_tax
                
    f = open(tax_genomad_tax_output, 'w')
    for gn in gn2lca_tax:
        f.write(f'{gn}\t{gn2lca_tax[gn]}\n')
    f.close()  
   
genomad_summary_file, vRhyme_best_bin_dir_modified, vRhyme_unbinned_viral_gn_dir, tax_genomad_tax_output = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
get_tax_from_genomad_result(genomad_summary_file, vRhyme_best_bin_dir_modified, vRhyme_unbinned_viral_gn_dir, tax_genomad_tax_output)