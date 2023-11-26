This folder contains pair-ended reads, metagenomic assembly, and custom MAGs for testing:

```bash
# TEST 1: Only use reads and metagenomic assembly
# example code for testing:
ViWrap run  --input_metagenome test_metaG.fasta \
            --input_reads reads_1.fastq.gz,reads_2.fastq.gz \
            --out_dir  test_metaG_ViWrap_out \
            --db_dir /storage1/data11/ViWrap/ViWrap_db \ # Change according to your case
            --identify_method vb-vs \
            --conda_env_dir /slowdata/yml_environments \ # Change according to your case
            --threads 20 \
            --input_length_limit 5000

# The total running time for TEST 1 is about 2 hrs  


# TEST 2: Use reads, metagenomic assembly, and custom MAGs (binned from the same metagenome)
# example code for testing:
ViWrap run  --input_metagenome test_metaG.fasta \
            --input_reads reads_1.fastq.gz,reads_2.fastq.gz \
            --out_dir  test_metaG_ViWrap_out \
            --db_dir /storage1/data11/ViWrap/ViWrap_db \ # Change according to your case
            --identify_method vb \
            --conda_env_dir /slowdata/yml_environments \ # Change according to your case
            --threads 20 \
            --input_length_limit 5000 \
            --custom_MAGs_dir custom_MAGs_dir \
            --iPHoP_db_custom iPHoP_db_custom
            
# The total running time for TEST 2 is about 19 hrs    
```

