#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    warnings.filterwarnings("ignore")
    from collections import defaultdict
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 
    
     
def get_tax_from_vcontact2_result(genome_by_genome_file, IMGVR_db_map, output):
    # Step 1 Get ref_gn2tax dict
    ref_gn2tax = {} # gn => tax
    with open(IMGVR_db_map, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            gn, tax = line.split(',')[1], line.split(',')[2]
            if gn not in ref_gn2tax:
                ref_gn2tax[gn] = tax
    
    # Step 2 Store VC2gns 
    VC2gns = defaultdict(list)  # VC => [gns]
    with open(genome_by_genome_file, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            line = line.rstrip('\n')
            split_line = line.split(',')
            VC, gn, genus_confidence_score = split_line[3], split_line[0], split_line[9] # gn here is actually the old name - the scaffold header
            if VC != '' and gn != 'contig_id' and genus_confidence_score != '' and float(genus_confidence_score) >= 0.9:
                VC2gns[VC].append(gn)                   
                        
    # Step 3 Get consensus tax
    bin2consensus_tax = {} # bin => consensus_tax
    for VC in VC2gns:
        gns = VC2gns[VC]
        consensus_tax = '' # Consensus tax for this VC (genus level), the number of taxes should be 1
        tax_set = set() # Store all unique tax
        for gn in gns:
            if gn in ref_gn2tax:
                tax = ref_gn2tax[gn]
                tax_set.add(tax)
        if len(tax_set) == 1:
            consensus_tax = list(tax_set)[0]
            
        if consensus_tax != '':
            for gn in gns:
                if gn not in ref_gn2tax:
                    bin2consensus_tax[gn] = consensus_tax
                
    # Step 4 Write to output
    f = open(output, 'w')
    for bin_name in bin2consensus_tax:
        f.write(f'{bin_name}\t{bin2consensus_tax[bin_name]}\n')
    f.close()     
   
   
genome_by_genome_file, IMGVR_db_map, output = sys.argv[1], sys.argv[2], sys.argv[3]
get_tax_from_vcontact2_result(genome_by_genome_file, IMGVR_db_map, output)    