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
    Realm = Kingdom = Phylum = Class = Order = Family = Genus = Species = {}
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
        lca_list = [Realm.keys[0]]
        
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) != 1:
        lca_list = [Realm.keys[0], Kingdom.keys[0]]
            
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) != 1:
        lca_list = [Realm.keys[0], Kingdom.keys[0], Phylum.keys[0]]

    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) != 1:
        lca_list = [Realm.keys[0], Kingdom.keys[0], Phylum.keys[0], Class.keys[0]]
            
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) == 1 and len(Family) != 1 :
        lca_list = [Realm.keys[0], Kingdom.keys[0], Phylum.keys[0], Class.keys[0], Order.keys[0]]            
            
    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) == 1 and len(Family) == 1 and len(Genus) != 1:
        lca_list = [Realm.keys[0], Kingdom.keys[0], Phylum.keys[0], Class.keys[0], Order.keys[0], Family.keys[0]]

    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) == 1 and len(Family) == 1 and len(Genus) == 1 and len(Species) != 1:
        lca_list = [Realm.keys[0], Kingdom.keys[0], Phylum.keys[0], Class.keys[0], Order.keys[0], Family.keys[0], Genus.keys[0]]

    if len(Realm) == 1 and len(Kingdom) == 1 and len(Phylum) == 1 and len(Class) == 1 and len(Order) == 1 and len(Family) == 1 and len(Genus) == 1 and len(Species) == 1:
        lca_list = [Realm.keys[0], Kingdom.keys[0], Phylum.keys[0], Class.keys[0], Order.keys[0], Family.keys[0], Genus.keys[0], Species.keys[0]]
    
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
                
def integrate_all_taxonomical_results(viwrap_outdir, genus_cluster_info, tax_classification_result):        
    tax_refseq_output = f'{viwrap_outdir}/tax_refseq_output.txt'
    tax_vog_output = f'{viwrap_outdir}/tax_vog_output.txt'
    tax_vcontact2_output = f'{viwrap_outdir}/tax_vcontact2_output.txt'
    
    tax_refseq = store_tax_output(tax_refseq_output)
    tax_vog = store_tax_output(tax_vog_output)
    tax_vcontact2 = store_tax_output(tax_vcontact2_output)
    
    gn2tax = {} # gn => [tax, method]
    # Step 1 Store taxonomic classification result from three methods
    for gn in tax_refseq:
        tax =  tax_refseq[gn]
        method = 'NCBI RefSeq viral protein searching'
        taxs = tax.split(';')
        for i in range(len(taxs)):
            if taxs[i] == '':
                taxs[i] = 'NA'
                
        tax = ';'.join(taxs)
        gn2tax[gn] = [tax, method]
        
    for gn in tax_vog:    
        tax = tax_vog[gn]
        method = 'marker VOG HMM searching'
        if gn not in gn2tax: # NCBI RefSeq viral protein searching method has the higher priority
            # Add genus and species into the tax (of course both are 'NA')
            taxs = tax.split(';')
            taxs.append('NA')
            taxs.append('NA')
            tax = ';'.join(taxs)
            gn2tax[gn] = [tax, method]  

    for gn in tax_vcontact2:    
        tax = tax_vcontact2[gn]
        method = 'vContact2 clustering'
        if gn not in gn2tax: # vContact2 clustering method has the lowest priority
            # Add genus and species into the tax (of course both are 'NA')
            gn2tax[gn] = [tax, method]  
            
                
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


viwrap_outdir, genus_cluster_info, tax_classification_result = sys.argv[1], sys.argv[2], sys.argv[3]
integrate_all_taxonomical_results(viwrap_outdir, genus_cluster_info, tax_classification_result)    