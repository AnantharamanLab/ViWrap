#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    import pysam
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call    
    warnings.filterwarnings("ignore")
    from pathlib import Path
    import pyfastx # For fastq and fasta reading and parsing 
    import pandas as pd
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)

def make_bowtie2_idx(fasta, working_dir, num_threads):
    file_name = Path(fasta).stem
    index_name = file_name + ".bowtie2_idx"
    
    fa = pyfastx.Fasta(fasta)
    fasta_size = fa.size
    
    if fasta_size <= 4000000000:
        # Indexing the reference sequence 
        indexing_cmd = f'bowtie2-build {fasta} {working_dir}/{index_name} --threads {num_threads} --quiet 1> /dev/null'
        os.system(indexing_cmd)
    else:
        # Indexing the reference sequence 
        indexing_cmd = f'bowtie2-build --large-index {fasta} {working_dir}/{index_name} --threads {num_threads} --quiet 1> /dev/null'
        os.system(indexing_cmd)
    
def run_bowtie2(fasta, input_read_pair, working_dir, sam_name, num_threads):
    file_name = Path(fasta).stem
    index_name = file_name + ".bowtie2_idx"
    
    # Mapping 
    mapping_cmd = f'bowtie2 -x {working_dir}/{index_name} -1 {input_read_pair.split(",")[0]} -2 {input_read_pair.split(",")[1]} -S {working_dir}/{sam_name}.sam -p {int(num_threads)} --no-unal --quiet --mm 1> /dev/null'
    os.system(mapping_cmd) 
    
def run_minimap2(fasta, input_reads, working_dir, sam_name, input_reads_type, num_threads):
    input_reads_type_map = {'pacbio':'map-pb', 'pacbio_hifi':'map-hifi', 'pacbio_asm20':'asm20', 'nanopore':'map-ont'}
    ax_input = input_reads_type_map[input_reads_type]
    
    # Mapping
    mapping_cmd = f'minimap2 -ax {ax_input} {fasta} {input_reads} -t {int(num_threads)} > {working_dir}/{sam_name}.sam 2> /dev/null' 
    os.system(mapping_cmd)  
    
def run_consent(input_reads, input_reads_type, num_threads):
    num_threads = int(num_threads)
    input_reads_type_map = {'pacbio':'PB', 'pacbio_hifi':'PB', 'pacbio_asm20':'PB', 'nanopore':'ONT'}
    reads_type = input_reads_type_map[input_reads_type]
    
    out_fasta_file = ''
    if '.gz' not in input_reads:
        out_fasta_file = input_reads.replace('.fastq', '.corrected.fasta', 1)
    elif '.gz' in input_reads:
        out_fasta_file = input_reads.replace('.fastq.gz', '.corrected.fasta', 1)
    
    # Correcting
    correcting_cmd = f'CONSENT-correct --in {input_reads} --out {out_fasta_file} --type {reads_type} -j {num_threads} 1> /dev/null'
    os.system(correcting_cmd)

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
        
def convert_sam_to_sorted_bam(input_sam_file, num_threads):
    # Open the SAM file in reading mode
    samfile = pysam.AlignmentFile(input_sam_file, "r")

    # Create a BAM file in writing mode
    out_bam_file = input_sam_file.replace('.sam', '.bam', 1)
    bamfile = pysam.AlignmentFile(out_bam_file, "wb", template=samfile)

    # Iterate through the records in the SAM file and write them to the BAM file
    for record in samfile:
        bamfile.write(record)

    # Close the files
    samfile.close()
    bamfile.close()

    # Sort the BAM file
    out_sorted_bam_file = input_sam_file.replace('.sam', '.sorted.bam', 1)
    pysam.sort("-o", out_sorted_bam_file, out_bam_file) 

def filter_sorted_bam(out_sorted_bam_file, filtered_bam_file, reads_mapping_identity_cutoff, aligned_length, threads):   
    reads_mapping_identity_cutoff = int(float(reads_mapping_identity_cutoff) * 100)
    threads = int(threads)
    filter_cmd = f'coverm filter --bam-files {out_sorted_bam_file} --output-bam-files {filtered_bam_file} --min-read-aligned-length {aligned_length} --min-read-percent-identity {reads_mapping_identity_cutoff} --threads {threads}'
    os.system(filter_cmd)
        
def mapping_metaG_reads(viral_scaffold, metagenomic_scaffold, metaG_reads, mapping_result_dir, input_reads_type, reads_mapping_identity_cutoff, threads):
    threads = int(threads)
    if input_reads_type == 'illumina':
        # Step 1 Run Bowtie2
        os.mkdir(mapping_result_dir)
        metaG_reads_list = metaG_reads.split(',')
        sam_names = []
        if len(metaG_reads_list) / 2 == 1:
            sam_name = Path(metaG_reads_list[0]).stem.rsplit('_', 1)[0]
            sam_names.append(sam_name)
            each_metaG_read_pair = ','.join(metaG_reads_list)
            make_bowtie2_idx(metagenomic_scaffold, mapping_result_dir, threads)
            run_bowtie2(metagenomic_scaffold, each_metaG_read_pair, mapping_result_dir, sam_name, threads)
        elif len(metaG_reads_list) / 2 >= 2 and len(metaG_reads_list) % 2 == 0:
            make_bowtie2_idx(metagenomic_scaffold, mapping_result_dir, threads)
            for i in range(0, len(metaG_reads_list), 2):
                j = i + 1
                sam_name = Path(metaG_reads_list[i]).stem.rsplit('_', 1)[0]
                sam_names.append(sam_name)
                each_metaG_read_pair = f'{metaG_reads_list[i]},{metaG_reads_list[j]}'
                run_bowtie2(metagenomic_scaffold, each_metaG_read_pair, mapping_result_dir, sam_name, threads)
        else:
            sys.exit('You input reads are not in pairs, please check') 
        
        # Step 2 Filter sam file
        for sam_name in sam_names:
            input_sam_file = f'{mapping_result_dir}/{sam_name}.sam'
            convert_sam_to_sorted_bam(input_sam_file, threads)
            out_bam_file = input_sam_file.replace('.sam', '.bam', 1)
            out_sorted_bam_file = input_sam_file.replace('.sam', '.sorted.bam', 1)
            filtered_bam_file = input_sam_file.replace('.sam', '.filtered.bam', 1)
            aligned_length = 50
            filter_sorted_bam(out_sorted_bam_file, filtered_bam_file, reads_mapping_identity_cutoff, aligned_length, threads)
            os.system(f'rm {input_sam_file} {out_bam_file} {out_sorted_bam_file}')    
        
        # Step 3 Get coverage
        bam_files = ''
        bam_files_list = []
        for sam_name in sam_names:
            bam_files_list.append(f'{mapping_result_dir}/{sam_name}.filtered.bam')
        bam_files = ' '.join(bam_files_list)
        
        os.system(f'coverm contig --methods metabat --bam-files {bam_files} --threads {threads} > {mapping_result_dir}/all_coverm_raw_result.txt')
        
        # Step 4 Parse all_coverm_raw_result.txt
        coverm_raw_table = pd.read_csv(f'{mapping_result_dir}/all_coverm_raw_result.txt', sep = '\t')
        coverm_raw_table_subset = coverm_raw_table.drop(['contigLen', 'totalAvgDepth'], axis = 1)
        
        dict_virus_rename = {} # old_name => new_name
        viral_seq = store_seq(viral_scaffold)
        for header in viral_seq:
            new_name = header.replace('>', '', 1)
            old_name = ''
            if '||' in new_name:
                old_name = new_name.rsplit('||', 1)[0]
            elif '_fragment_' in new_name:
                old_name = new_name.rsplit('_fragment_', 1)[0]
            else:
                old_name = new_name
            dict_virus_rename[old_name] = new_name   
        
        coverm_raw_table_subset.replace({"contigName": dict_virus_rename},inplace = True)
        coverm_raw_table_subset.to_csv(f'{mapping_result_dir}/vRhyme_input_coverage.txt', sep='\t', index=False)
    elif input_reads_type == 'pacbio' or input_reads_type == 'pacbio_hifi' or input_reads_type == 'pacbio_asm20' or input_reads_type == 'nanopore':
        # Step 1 Run minimap2
        os.mkdir(mapping_result_dir)      
        metaG_reads_list = []
        if ',' in metaG_reads:
            metaG_reads_list = metaG_reads.split(',')
        else:
            metaG_reads_list = [metaG_reads]
        
        sam_names = []
        for each_read in metaG_reads_list:
            sam_name = Path(each_read).stem
            sam_names.append(sam_name)
            run_consent(each_read, input_reads_type, threads)
            corrected_fasta_file = ''
            if '.gz' not in each_read:
                corrected_fasta_file = each_read.replace('.fastq', '.corrected.fasta', 1)
            elif '.gz' in each_read:
                corrected_fasta_file = each_read.replace('.fastq.gz', '.corrected.fasta', 1)
            run_minimap2(metagenomic_scaffold, corrected_fasta_file, mapping_result_dir, sam_name, input_reads_type, threads)

        # Step 2 Filter sam file
        for sam_name in sam_names:
            input_sam_file = f'{mapping_result_dir}/{sam_name}.sam'
            convert_sam_to_sorted_bam(input_sam_file, threads)
            out_bam_file = input_sam_file.replace('.sam', '.bam', 1)
            out_sorted_bam_file = input_sam_file.replace('.sam', '.sorted.bam', 1)
            filtered_bam_file = input_sam_file.replace('.sam', '.filtered.bam', 1)
            aligned_length = 500
            filter_sorted_bam(out_sorted_bam_file, filtered_bam_file, reads_mapping_identity_cutoff, aligned_length, threads)
            os.system(f'rm {input_sam_file} {out_bam_file} {out_sorted_bam_file}')   
        
        # Step 3 Get coverage
        bam_files = ''
        bam_files_list = []
        for sam_name in sam_names:
            bam_files_list.append(f'{mapping_result_dir}/{sam_name}.filtered.bam')
        bam_files = ' '.join(bam_files_list)
        
        os.system(f'coverm contig --methods metabat --bam-files {bam_files} --threads {threads} > {mapping_result_dir}/all_coverm_raw_result.txt')   

        # Step 4 Parse all_coverm_raw_result.txt
        coverm_raw_table = pd.read_csv(f'{mapping_result_dir}/all_coverm_raw_result.txt', sep = '\t')
        coverm_raw_table_subset = coverm_raw_table.drop(['contigLen', 'totalAvgDepth'], axis = 1)
        
        dict_virus_rename = {} # old_name => new_name
        viral_seq = store_seq(viral_scaffold)
        for header in viral_seq:
            new_name = header.replace('>', '', 1)
            old_name = ''
            if '||' in new_name:
                old_name = new_name.rsplit('||', 1)[0]
            elif '_fragment_' in new_name:
                old_name = new_name.rsplit('_fragment_', 1)[0]
            else:
                old_name = new_name
            dict_virus_rename[old_name] = new_name   
        
        coverm_raw_table_subset.replace({"contigName": dict_virus_rename},inplace = True)
        coverm_raw_table_subset.to_csv(f'{mapping_result_dir}/vRhyme_input_coverage.txt', sep='\t', index=False)        
    
    
viral_scaffold, metagenomic_scaffold, metaG_reads, mapping_result_dir, input_reads_type, reads_mapping_identity_cutoff, threads = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]
mapping_metaG_reads(viral_scaffold, metagenomic_scaffold, metaG_reads, mapping_result_dir, input_reads_type, reads_mapping_identity_cutoff, threads)       
    
    
    