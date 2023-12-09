#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from genomad import database, mmseqs2, utils
    from genomad._paths import GenomadOutputs    
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1) 
    
def run_genomad(metagenomic_scaffold, viwrap_outdir, threads, db_dir):
    cmd = f'genomad end-to-end {metagenomic_scaffold}  {viwrap_outdir}/genomad_output {db_dir}/genomad_db -t {threads} -q --cleanup --conservative-taxonomy 1> /dev/null'
    os.system(cmd)     
    
def identify_integrases(input_path, output_path, database_path, threads, sensitivity, evalue):
    """
    This function identifies integrases in genomic sequences using the MMseqs2 tool and the geNomad database.

    Args:
    - input_path (Path): Path to the annotated proteins output file.
    - output_path (Path): Path to the directory where the MMseqs2 outputs will be saved.
    - database_path (Path): Path to the geNomad database.
    - threads (int): Number of threads to use for MMseqs2 execution.
    - sensitivity (float): Sensitivity setting for MMseqs2.
    - evalue (float): E-value threshold for MMseqs2.

    Returns:
    None. The function writes outputs to files in the specified directory.
    """

    # Define the outputs object for managing output file paths
    outputs = GenomadOutputs(input_path.stem, output_path)

    # Initialize MMseqs2 with the geNomad database for integrase identification
    mmseqs2_obj = mmseqs2.MMseqs2(outputs.find_proviruses_mmseqs2_output, outputs.find_proviruses_mmseqs2_dir, input_path, database.Database(database_path), use_integrase_db=True)

    # Ensure the MMseqs2 directory and its parents are created
    if not outputs.find_proviruses_mmseqs2_dir.exists():
        outputs.find_proviruses_mmseqs2_dir.mkdir(parents=True, exist_ok=True)

    # Run MMseqs2 for integrase identification
    mmseqs2_obj.run_mmseqs2(threads, sensitivity, evalue, 0)    
    
    # Change the output dir and file names
    new_output_file = f"{output_path}/{input_path.stem}_find_proviruses/{input_path.stem}_provirus_mmseqs2.tsv"
    mv_cmd = f"mv {new_output_file} {output_path}/{input_path.stem}_integrase_mmseqs2.tsv"
    os.system(mv_cmd)
    delete_old_output_cmd = f"rm -rf {output_path}/{input_path.stem}_find_proviruses/"
    os.system(delete_old_output_cmd)   
    
# Take the input arguments
metagenomic_scaffold, viwrap_outdir, threads, db_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]  

# Run geNomad
run_genomad(metagenomic_scaffold, viwrap_outdir, threads, db_dir)    

# Call the function to identify integrases
input_path = f"{viwrap_outdir}/genomad_output/{Path(metagenomic_scaffold).stem}_summary/{Path(metagenomic_scaffold).stem}_virus_proteins.faa"
output_path = f"{viwrap_outdir}/genomad_output/integrase_result_dir"
database_path = f"{db_dir}/genomad_db"
sensitivity = 8.2  # Sensitivity for MMseqs2
evalue = 0.001  # E-value threshold for MMseqs2
identify_integrases(input_path, output_path, database_path, threads, sensitivity, evalue)