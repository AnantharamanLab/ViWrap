#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings("ignore")   
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 

def add_custom_MAGs_to_host_db__make_gtdbtk_results(viwrap_outdir, custom_MAGs_dir, threads):
    fasta_files = glob(os.path.join(custom_MAGs_dir, '*.fasta')) 
    if not fasta_files:
        # Step 1 Check the extension of input genomes
        sys.exit(f'Please make sure your there are input MAGs in {custom_MAGs_dir} and all of them end with ".fasta"')
    else:
        os.system(f'gtdbtk de_novo_wf --genome_dir {custom_MAGs_dir} --bacteria --outgroup_taxon p__Patescibacteria --out_dir {viwrap_outdir}/07_iPHoP_outdir/custom_MAGs_GTDB-tk_results --cpus {threads} --force --extension fasta 1> /dev/null')
        os.system(f'gtdbtk de_novo_wf --genome_dir {custom_MAGs_dir} --archaea --outgroup_taxon p__Altarchaeota --out_dir {viwrap_outdir}/07_iPHoP_outdir/custom_MAGs_GTDB-tk_results --cpus {threads} --force --extension fasta 1> /dev/null')
                
viwrap_outdir, custom_MAGs_dir, threads = sys.argv[1], sys.argv[2], sys.argv[3]
add_custom_MAGs_to_host_db__make_gtdbtk_results(viwrap_outdir, custom_MAGs_dir, threads)        