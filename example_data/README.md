This folder contains input metagenomic assembly and pair-ended reads for testing:

```bash
# Example code for testing:
ViWrap run  --input_metagenome test_metaG.fasta \
            --input_reads reads_1.fastq.gz,reads_2.fastq.gz \
            --out_dir  test_metaG_ViWrap_out \
            --db_dir /storage1/data11/ViWrap/ViWrap_db \ # Change according to your case
            --identify_method vb-vs \
            --conda_env_dir /slowdata/yml_environments \ # Change according to your case
            --threads 10 \
            --input_length_limit 5000
```

