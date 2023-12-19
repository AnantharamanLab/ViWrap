#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 
    
   
def get_LCA_from_taxs(taxs):
    Realm = {}
    Kingdom = {} 
    Phylum = {} 
    Class = {} 
    Order = {} 
    Family = {} 
    Genus = {} 
    Species = {}
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
            
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) == 1 and len(Family) != 1 :
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

def store_tax_output(input_file):
    tax_output = {} # gn => tax
    with open(input_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            gn, tax = tmp[0], tmp[1]
            tax_output[gn] = tax
    return tax_output  

def get_lowest_non_na_rank(tax):
    # Split the taxonomy string and reverse it to start checking from the lowest rank
    tax_ranks = tax.split(';')[::-1]
    for i, rank in enumerate(tax_ranks):
        if rank != 'NA':
            return i  # Return the index of the lowest non-NA rank
    return float('inf')  # Return infinity if all ranks are 'NA'    
                
def integrate_all_taxonomical_results(identify_method, viwrap_outdir, genus_cluster_info, tax_classification_result):        
    tax_refseq_output = f'{viwrap_outdir}/tax_refseq_output.txt'
    tax_vog_output = f'{viwrap_outdir}/tax_vog_output.txt'
    tax_vcontact2_output = f'{viwrap_outdir}/tax_vcontact2_output.txt'
    tax_genomad_output = f'{viwrap_outdir}/tax_genomad_tax_output.txt'
    
    tax_refseq = store_tax_output(tax_refseq_output)
    tax_vog = store_tax_output(tax_vog_output)
    tax_vcontact2 = store_tax_output(tax_vcontact2_output)
    tax_genomad = {}
    if identify_method == 'genomad':
        tax_genomad = store_tax_output(tax_genomad_output)
        
    if len(tax_refseq) > 0: # Modify the tax of tax_refseq 
        for gn in tax_refseq:
            tax = tax_refseq[gn]
            taxs = tax.split(';')
            for i in range(len(taxs)):
                if taxs[i] == '':
                    taxs[i] = 'NA'                
            tax = ';'.join(taxs)
            tax_refseq[gn] = tax

    if len(tax_vog) > 0: # Modify the tax of tax_vog     
        for gn in tax_vog:    
            tax = tax_vog[gn]
            taxs = tax.split(';')
            taxs.append('NA')
            taxs.append('NA')
            tax = ';'.join(taxs)
            tax_vog[gn] = tax

    # Step 1 Get gn2tax dict    
    ## Initialize the dictionary if it's not already present
    gn2tax = {}

    ## Process each genome name for each method and update the dictionary based on the new priority logic
    all_gns = set(tax_refseq.keys()) | set(tax_vog.keys()) | set(tax_vcontact2.keys()) | set(tax_genomad.keys())
    for gn in all_gns:
        best_method = None
        best_rank = float('inf')

        ### Check each method
        for method, tax_dict in [('NCBI RefSeq viral protein searching', tax_refseq),
                                 ('marker VOG HMM searching', tax_vog),
                                 ('vContact2 clustering', tax_vcontact2),
                                 ('geNomad taxonomy', tax_genomad)]:
            tax = tax_dict.get(gn, 'NA;NA;NA;NA;NA;NA;NA;NA')  # Default to 'NA' for all ranks if not found
            rank = get_lowest_non_na_rank(tax)
            if rank < best_rank:
                best_rank = rank
                best_method = method
                best_tax = tax

        ### Update the gn2tax dictionary with the best method
        if best_method:
            gn2tax[gn] = [best_tax, best_method]
          
    # Step 2 Store genus info
    genus2gns = {} # genus => [gns]
    with open(genus_cluster_info ,'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#'):
                genus = line.split(',', 1)[0]
                gns = line.split(',', 1)[1].split(';')
    
    # Step 3 Get taxonomy based on other members' taxonomy from each genus
    # Get into each genus cluster to see if any genomes have already got hits 
    # (only the hits of 'NCBI RefSeq viral protein searching' will be counted), 
    # then expand the tax to all the members within this genus cluster
    gn2tax_by_genus = {} # gn => [tax, method] 
    for genus in genus2gns:
        gns = genus2gns[genus]
        
        taxs = [] # The taxonomy collection for all the genomes within this genus
        for gn in gns:
            if gn in gn2tax:
                if gn2tax[gn][1] == 'NCBI RefSeq viral protein searching':
                    taxs.append(gn2tax[gn][0])
                    
        lca = 'NA;NA;NA;NA;NA;NA;NA;NA'; # The LCA for all taxonomic hits within this genus (should be above genus level)
        if taxs:
            lca = get_LCA_from_taxs(taxs)
            if lca != 'NA;NA;NA;NA;NA;NA;NA;NA':
                for gn in gns:
                    if gn not in gn2tax:
                        gn2tax_by_genus[gn] = [lca, 'Genus LCA assigning']
                        
    # Step 4 Write down final taxonomic classification result                   
    f = open(tax_classification_result, 'w')
    for gn in gn2tax:
        f.write(f'{gn}\t{gn2tax[gn][0]}\t{gn2tax[gn][1]}\n')
    for gn in gn2tax_by_genus:
        f.write(f'{gn}\t{gn2tax_by_genus[gn][0]}\t{gn2tax_by_genus[gn][1]}\n')
    f.close()   


identify_method, viwrap_outdir, genus_cluster_info, tax_classification_result = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
integrate_all_taxonomical_results(identify_method, viwrap_outdir, genus_cluster_info, tax_classification_result)    