This folder contains input metagenomic assembly and pair-ended reads for testing:



`ViWrap run --inp`ut_metagenome test_metaG.fasta \

​                   --input_reads reads_1.fastq.gz,reads_2.fastq.gz \
​                   `--out_dir  test_metaG_ViWrap_out\
​                   `

​                    --db_dir /storage1/data11/ViWrap/ViWrap_db \
​                   `--identify_method vb-vs \`
​                   `--conda_env_dir /slowdata/yml_environments \`
​                   `--threads 10 \`
​                 --input_length_limit 2000`