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
    mapping_cmd = f'bowtie2 -x {working_dir}/{index_name} -1 {input_read_pair.split(",")[0]} -2 {input_read_pair.split(",")[1]} -S {working_dir}/{sam_name}.sam -p {num_threads} --no-unal --quiet --mm 1> /dev/null'
    os.system(mapping_cmd)    
        
class FilterCoverage: 
# This chunk of scripts were copied from 'filter_coverage_file.py' (within our custom python scripts collection) 
# by Kristopher Kieft, 2022 at University of Wisconsin-Madison
    def __init__(self, sam, bam, out, percent, edit, threads, length, method, interest):
        self.sam = sam # input SAM file
        self.bam = bam # input BAM file (sorted or unsored)
        self.out = out # output filtered file
        self.interest = interest # file containing list of scaffolds of interest
        self.percent = percent # percent identity cutoff [0.97] (default)
        self.edit = edit # maximum edit distance (mismatch+gap+insert+delete)
        self.threads = threads # threads (SAM->BAM, BAM sorting) [1]
        self.length = length # minimum length per read [50]
        self.method = method # method of filtering, choose within percent identity cutoff and maximum edit distance
        if not self.out:
            self.set_out()
        
        self.percent = round(1.0 - self.percent,2)
        
        if self.interest:
            self.get_interest()
        
        if self.sam:
            self.aligned_sam()
        else:
            self.aligned_bam()
        
        self.sort_check()

        if self.method == 'p':
            self.main_percent()
        else:
            self.main_edist()
        
        if self.indexed:
            subprocess.run(f'rm {self.alignment}.bai', shell=True)

    
    def set_out(self):
        if self.method == 'p':
            base = f'pid{int(self.percent*100)}'
        else:
            base = f'dist{self.edit}'
        if self.sam:
            self.out = self.sam.rsplit('.',1)[0] + f'.{base}.sorted.bam'
        elif self.bam:
            self.out = self.bam.rsplit('.',1)[0] + f'.{base}.sorted.bam'
        if os.path.exists(self.out):
            print('\nError: output filtered BAM already exists. Exiting.')
            print(f'{self.out}\n')
            exit()
    
    def get_interest(self):
        with open(self.interest) as f:
            self.names = set(f.read().split('\n'))
    
    def aligned_sam(self):
        check_aligned = subprocess.check_output(f"samtools view -@ {self.threads} {self.sam} | head -n 1", shell=True)
        if len(check_aligned) == 0:
            print(f'\nNo aligned reads in {self.sam}.\n')
            exit()
        try:
            temp = self.sam.rsplit('/',1)[1]
            self.bam = temp.rsplit('.',1)[0] + '.bam'
        except IndexError:
            self.bam = self.sam.rsplit('.',1)[0] + '.bam'
        if os.path.exists(self.bam):
            print('\nError: looks like a BAM file already exists. Exiting.')
            print(f'{self.bam}\n')
            exit()
        subprocess.run(f"samtools view -@ {self.threads} -S -b {self.sam} > {self.bam}", shell=True)
    
    def aligned_bam(self):
        check_aligned = subprocess.check_output(f"samtools view -@ {self.threads} {self.bam} | head -n 1", shell=True)
        if len(check_aligned) == 0:
            print(f'\nNo aligned reads in {self.bam}.\n')
            exit()

    def sort_check(self):
        sort_check = False
        try:
            check = subprocess.check_output(f"samtools view -@ {self.threads} -H {self.bam} | grep '@HD'", shell=True)
            if "coordinate" in str(check):
                sort_check = True
        except Exception:
            # no @HD line
            pass
        if not sort_check:
            sbam = self.bam.rsplit('.',1)[0] + '.sorted.bam'
            if os.path.exists(sbam):
                sbam = self.bam.rsplit('.',1)[0] + '.py-sorted.bam'

            subprocess.run(f"samtools sort -@ {self.threads} -o {sbam} {self.bam}", shell=True)
            self.alignment = sbam
        else:
            self.alignment = self.bam

        self.indexed = False
        if not os.path.exists(f'{self.alignment}.bai'):
            self.indexed = True
            subprocess.run(f'samtools index {self.alignment}', shell=True)

    def align_id(self, ed, rl, x):
        if ed/rl <= self.percent and rl >= self.length:
            self.outfile.write(x)

    def align_e(self, ed, rl, x):
        if ed <= self.edit and rl >= self.length:
            self.outfile.write(x)

    def main_percent(self):

        bamfile = pysam.AlignmentFile(self.alignment, "rb")
        self.outfile = pysam.AlignmentFile(self.out, "wb", template=bamfile)

        if self.interest:
            for x in bamfile.fetch(until_eof=True):
                genome = x.reference_name
                if genome in self.names:
                    rl = x.query_length
                    ed = 0
                    for t in x.tags:
                        if t[0] == 'NM':
                            ed = t[1]
                            self.align_id(ed, rl, x)
                            break

        else:
            for x in bamfile.fetch(until_eof=True):
                rl = x.query_length
                ed = 0
                for t in x.tags:
                    if t[0] == 'NM':
                        ed = t[1]
                        self.align_id(ed, rl, x)
                        break
    
        bamfile.close()
        self.outfile.close()
    
    def main_edist(self):

        bamfile = pysam.AlignmentFile(self.alignment, "rb")
        self.outfile = pysam.AlignmentFile(self.out, "wb", template=self.alignment)

        if self.interest:
            for x in bamfile.fetch(until_eof=True):
                genome = x.reference_name
                if genome in self.names:
                    rl = x.query_length
                    ed = False
                    for t in x.tags:
                        if t[0] == 'NM':
                            ed = t[1]
                            self.align_e(ed, rl, x)
                            break

        else:
            for x in bamfile.fetch(until_eof=True):
                rl = x.query_length
                ed = False
                for t in x.tags:
                    if t[0] == 'NM':
                        ed = t[1]
                        self.align_e(ed, rl, x)
                        break
    

        bamfile.close()
        self.outfile.close()        
        
def mapping_metaG_reads(metagenomic_scaffold, metaG_reads, mapping_result_dir, threads):
    threads = int(threads)
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
        FilterCoverage(f'{mapping_result_dir}/{sam_name}.sam', '', f'{mapping_result_dir}/{sam_name}.filtered.bam', 0.97, '', threads, 50, 'p', '')
    # sam, bam, out, percent, edit, threads, length, method, interest
    
    # sam: input SAM file
    # bam: input BAM file (sorted or unsored)
    # out: output filtered file
    # percent: percent identity cutoff [default: 0.97]
    # edit: maximum edit distance (mismatch+gap+insert+delete)
    # threads: threads (SAM->BAM, BAM sorting)
    # length: minimum length per read [default: 50]
    # method: method of filtering, choose within percent identity cutoff ('p') and maximum edit distance ('e')
    # interest: file containing list of scaffolds of interest
    
    os.system(f'rm {metagenomic_scaffold}.fxi')
    for sam_name in sam_names:
        os.system(f'rm {sam_name}.bam {sam_name}.sorted.bam {mapping_result_dir}/{sam_name}.sam')    
    
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
    coverm_raw_table_subset.to_csv(f'{mapping_result_dir}/vRhyme_input_coverage.txt', sep='\t', index=False)
    
    
metagenomic_scaffold, metaG_reads, mapping_result_dir, threads = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
mapping_metaG_reads(metagenomic_scaffold, metaG_reads, mapping_result_dir, threads)       
    
    
    