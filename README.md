<p align="center">
  <img width="350"  src="https://github.com/AnantharamanLab/ViWrap/blob/main/images/ViWrap2.jpg">
</p>

# **ViWrap**
### ViWrap: A modular pipeline to identify, bin, classify, and predict viral-host relationships for viruses from metagenomes 

```
Oct 2022   
Zhichao Zhou  
zzhou388@wisc.edu & zczhou2017@gmail.com  
Anantharaman Lab
Department of Bacteriology
University of Wisconsin-Madison  
```

## Current Version
ViWrap v1.2.1 

<img src="https://github.com/AnantharamanLab/ViWrap/blob/main/images/ViWrap-Flowchart.jpg">

## Citation
If you find ViWrap useful please consider citing our manuscript on XXX (still in draft):  

```tex
Zhou, Z et al. ViWrap: A modular pipeline to identify, bin, classify, and predict viral-host relationships for viruses from metagenomes. XXX 2023
```

______
## Table of Contents:
1. [Updates](#updates)
2. [Program Description](#program)
3. [Installation](#install)
4. [Settings](#settings)
5. [Running ViWrap](#run)
    * ViWrap tasks
    * Flag explanations
6. [Output Explanations](#out)
7. [Contact](#contact)

______
## Updates for v1.2.1 (Jan 2023): <a name="updates"></a>

* Provide "AMG_results" in "ViWrap_summary_outdir" to include AMG statistics, AMG protein details, and AMG protein sequences.
* Check and update the vContact2 conda environment.
* Correct function "get_virus_genome_annotation_result" in "module.py" ("==" changed to "in" in outside elif)
* Correct thread number mistake in "run_CheckV.py" script
* Add long reads mapping function (for nanopore or pacbio reads)
* Solve VirSorter 2 conda environment and database issues (set numpy to v1.20)
* Update bam filtering steps
* Update GTDB-Tk to v2.1.1 and db to release 207 v2
* Change the output place of "iPHoP_db_custom"
* Update the "master_downloader.py" to add checkV DB setting up command line
* Update the "module.py" to change long headers into short headers in the output virus faa files (within the folders of "02_vRhyme_outdir/vRhyme_best_bins_fasta_modified" and "08_ViWrap_summary_outdir/Virus_genomes_files"). And correct the ffn files error.
* Update "VOG_marker_table.txt" to remove "Caudovirales" from the tax column
* Update the IMG/VR v3 to v4
* Update the KEGG metabolism calculation issue (denominator sometimes could be "0") in the function "generate_result_visualization_inputs" within the "module.py" script
* Update the vs virus identifying method - delete the 2nd vs and CheckV steps and add a pre-check step for input scaffold length limit when using vs method

## Updates for v1.2.0 (Oct 2022): <a name="updates"></a>

* Modify "parse_dRep" function (module.py) so that it can also process the result even though dRep running is not successful.
* Correct the faa protein ID issue in the "get_virus_genome_annotation_result" function (module.py).
* Add "-red 5" (the maximum number of redundant proteins per bin should be less than 5) to the "run_vRhyme.py" script.
* Correct AMG KO statistic mistake in "master_run.py"; now only AMG KOs will be included in the final "Virus_summary_info.txt" file.
* Add "run_wo_reads" option to allow directly identifying viral scaffolds from metagenomes/genomes without the usage of metagenomic/genomic reads.
* Introduce "vb-vs" and "vb-vs-dvf" to the "identify_method" option to allow using the overlapped virus identification results for downstream analysis.
* Provide flowchart.
* Provide "scf2lytic_or_lyso.summary.txt" for the VIBRANT result, "vRhyme_best_bin_lytic_and_lysogenic_info.txt" and 
  "vRhyme_best_bin_scaffold_complete_info.txt" for vRhyme generated vRhyme_best_bins. Then, based on the formulas of "Lytic and lysogenic viruses/scaffolds" (See "Notes - Lytic and lysogenic viruses can bin together" in "Output Explanations"), make modified vRhyme_best_bins (some bins are split into scaffolds). The corresponding downstream analysis has also been changed accordingly.
* Add visualization output directory



## Updates for v1.1.0 (Sep 2022): <a name="updates"></a>

* Add two additional bypasses for identifying viruses from metagenomes: vs - VirSorter2 and CheckV; dvf - DeepVirFinder besides the default pass of vb - VIBRANT.
* set up the VirSorter2, CheckV, and VIBRANT virus identification, curation, and annotation bypass mainly according to the SOP provided by Matthew Sullivan's group (OSU); the 'manual check' step changed to be 'automatic' with stringent criteria.
* Complete GitHub README.


## Updates for v1.0.0 (Sep 2022): <a name="updates"></a>
* Change the construction of scripts.
* Change the methods to feed arguments and run scripts in individual conda environments.

______
## Program Description <a name="program"></a>

#### **ViWrap Description**

ViWrap is a wrapper to identify, bin, classify, and predict host-viral relationships for viruses from metagenomes. It leverages the advantages of currently available virus analyzing tools and provides a quick, intuitive, one-step pipeline to get viral sequences and corresponding properties.

*Note*: ViWrap is built for studying viruses identified and reconstructed from metagenomes. It is suggested to be only used for viruses with hosts as bacteria, archaea, and microeukaryotes.

ViWrap is an integrated wrapper/pipeline, the main contributors of each virus identifying, binning, classifying, and viral host predicting software within it should be acknowledged (Citations and links are provided):

[VIBRANT](https://github.com/AnantharamanLab/VIBRANT): [link to online paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0)

```
Kieft, Kristopher, Zhichao Zhou, and Karthik Anantharaman. "VIBRANT: automated recovery, annotation and curation of microbial viruses, and evaluation of viral community function from genomic sequences." Microbiome 8, no. 1 (2020): 1-23. 
```

[VirSorter2](https://github.com/jiarong/VirSorter2): [link to online paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y)

```
Guo, Jiarong, Ben Bolduc, Ahmed A. Zayed, Arvind Varsani, Guillermo Dominguez-Huerta, Tom O. Delmont, Akbar Adjie Pratama et al. "VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses." Microbiome 9 (2021): 1-13.
```

[DeepVirFinder](https://github.com/jessieren/DeepVirFinder): [link to online paper](https://link.springer.com/article/10.1007/s40484-019-0187-4)

```
Ren, Jie, Kai Song, Chao Deng, Nathan A. Ahlgren, Jed A. Fuhrman, Yi Li, Xiaohui Xie, Ryan Poplin, and Fengzhu Sun. "Identifying viruses from metagenomic data using deep learning." Quantitative Biology 8 (2020): 64-77.
```

[vContact2](https://bitbucket.org/MAVERICLab/vcontact2/wiki/Home): [link to online paper](https://www.nature.com/articles/s41587-019-0100-8)

```
Bin Jang, Ho, Benjamin Bolduc, Olivier Zablocki, Jens H. Kuhn, Simon Roux, Evelien M. Adriaenssens, J. Rodney Brister et al. "Taxonomic assignment of uncultivated prokaryotic virus genomes is enabled by gene-sharing networks." Nature biotechnology 37, no. 6 (2019): 632-639.
```

[vRhyme](https://github.com/AnantharamanLab/vRhyme): [link to online paper](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac341/6584432)

```
Kieft, Kristopher, Alyssa Adams, Rauf Salamzade, Lindsay Kalan, and Karthik Anantharaman. "vRhyme enables binning of viral genomes from metagenomes." Nucleic Acids Research 50, no. 14 (2022): e83-e83.
```

[iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/) ([and software within it](https://bitbucket.org/srouxjgi/iphop/src/main/#markdown-header-citation)): [link to online paper](https://www.biorxiv.org/content/10.1101/2022.07.28.501908v1)

```
Roux, Simon, Antonio Pedro Camargo, Felipe Hernandes Coutinho, Shareef M. Dabdoub, Bas E. Dutilh, Stephen Nayfach, and Andrew Tritt. "iPHoP: an integrated machine-learning framework to maximize host prediction for metagenome-assembled virus genomes." bioRxiv (2022): 2022-07.
```




#### **ViWrap Features**
* Identify and annotate viruses by VIBRANT, VirSorter2, or DeepVirFinder
* Map reads (support both short and long reads) onto all input metagenome assemblies to get scaffold depth
* Bin vMAGs using vRhyme
* Classify vMAGs into genus using vContact2, into species using dRep
* Get virus quality using CheckV
* Get virus taxonomy using three approaches
* Get viral hosts using iPHoP 
* Get virus genome abundance
* Get all virus genome information

______
## Installation <a name="install"></a>

#### **Set up the conda environment for ViWrap**

Since ViWrap has many dependencies to be installed, it would be much easier to set up a conda environment instead of installing all dependencies in the global environment (make sure you have upfront conda installed on your server, i.e., [miniconda3](https://docs.conda.io/en/latest/miniconda.html) or anaconda3; we only suggest to run in version 3.0+ conda). There are 12 conda envs associated with ViWrap, so it may be useful to keep these associated together in a directory separate from your normal conda installation, which will be referred to as `/path/to/ViWrap_conda_environments`. **Note: ensure that wherever you install the ViWrap conda envs has at least 5-6 Gb of storage available.** ViWrap was tested extensively with an environment installed separately from the normal conda installation.

Choose one:
1. Install in separate directory:
   1. `conda create -c bioconda -c conda-forge -p /path/to/ViWrap_conda_environments/ViWrap python=3.8 biopython=1.80 mamba=1.3.0 numpy=1.24.2 pandas=1.5.3 pyfastx=0.8.4 matplotlib=3.6.3 seaborn=0.12.2 diamond=2.0.15 hmmer=3.3.2`
   2. `conda activate /path/to/ViWrap_conda_environments/ViWrap`

   Note: `/path/to/conda_environments` indicates the directory that you will need to use to store all conda environments for ViWrap
2. Install in normal conda folder
   1. `conda create -c bioconda -c conda-forge -n ViWrap python=3.8 biopython=1.80 mamba=1.3.0 numpy=1.24.2 pandas=1.5.3 pyfastx=0.8.4 matplotlib=3.6.3 seaborn=0.12.2 diamond=2.0.15 hmmer=3.3.2`
   2. `conda activate ViWrap`
   
   If you choose to proceed this route, you will just need to use the path to your ViWrap conda installation. It will look something like this: `/path/to/miniconda3/envs/ViWrap/`

#### GitHub installation

1. `git clone https://github.com/AnantharamanLab/ViWrap`
2. `cd ViWrap`
3. `chmod +x ViWrap scripts/*.py # Make all python scripts to be executable`
4. ``PATH=`pwd`:$PATH # Add ViWrap to the PATH, so it can be called elsewhere in a terminal``

#### Set up the other conda environments required by ViWrap

`ViWrap set_up_env --conda_env_dir /path/to/ViWrap_conda_environments`

This will take several minutes depending on your current internet speed. 

**Note**: `/path/to/ViWrap_conda_environments` can be set anywhere on your server to contain ViWrap conda environments. If you use `-n ViWrap` in the above step, you will input `/path/to/miniconda3/envs` to `--conda_env_dir` option.

ViWrap will use the "-p" or "--prefix" option to specify where to write the environment files:

(Details in https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#specifying-a-location-for-an-environment)

For example,`conda create --prefix /tmp/test-env python=3.8`
will create the environment named `/tmp/test-env` which resides in `/tmp/` instead of the default `.conda`.

The following 11 conda environments will be set up, the **estimated running time will be 20-30 minutes**, depending on your current internet speed:

```
121M	./ViWrap-VIBRANT
265M	./ViWrap-vContact2
100M	./ViWrap-dRep
75M	./ViWrap-Tax
205M	./ViWrap-iPHoP
1.9G	./ViWrap-vs2
272M	./ViWrap-vRhyme
76M	./ViWrap-CheckV
83M	./ViWrap-GTDBTk
155M	./ViWrap-Mapping
2.1G	./ViWrap-DVF
```

**Note:**  We have fixed the versions of Python modules to prevent potential errors caused by version upgrades. If you encounter any issues with these conda environments, please verify the module versions set in "scripts/master_set_up_env.py". 

#### Set up ViWrap database

```shell
ViWrap download --db_dir /path/to/ViWrap_db  --conda_env_dir /path/to/ViWrap_conda_environments
```

`/path/to/ViWrap_db` is the place you store the ViWrap database. Please make sure there is enough space to store the database (~280G at least). It will **take** **~3-4 hours** to set up well depending on your current internet speed. This is kind of tedious, however, you will only need to do this one time.

It contains the following 7 folders (call by `du -h --max-depth=1 ./` within the directory of "ViWrap_db"): 

```
11G	./VIBRANT_db
6.4G	./CheckV_db
114M	./DVF_db
581M	./Tax_classification_db
199G	./iPHoP_db
11G	./VirSorter2_db
49G	./GTDB_db
```

***Notes:*** 

1) Since some software (VirSorter2) needs to config the database address into the conda environment, it is suggest to first set up the environments, then set up the databases. 

2) Once you have replaced any conda environments, it is better to re-check/re-install the corresponding conda environments (especially for the case of VirSorter2)

3) Since some db folders are restricted within the creator's rights, if anyone else in the group who wants to use ViWrap, the db folder rights should be opened by using`chmod -R 777 ./`

#### See ViWrap help

1. `ViWrap -h`
2. `ViWrap run -h`

------

## Settings <a name="settings"></a>

#### Settings for v1.2.1

- **Conda env:**

  | Software   | Version | Software                        | Version     |
  | ---------- | ------- | ------------------------------- | ----------- |
  | VIBRANT    | v1.2.1  | dRep                            | v3.4.0      |
  | vRhyme     | v1.1.0  | DIAMOND (within ViWrap-Tax)     | v2.0.15     |
  | vContact2  | v0.11.0 | iPHoP                           | v1.2.0      |
  | CheckV     | v1.0.1  | GTDB-Tk                         | v2.1.1      |
  | VirSorter2 | v2.2.3  | Bowtie2 (within ViWrap-Mapping) | v2.4.5      |
  | Minimap2   | v2.24   | DeepVirFinder                   | v2020.11.21 |

- **db:**

  | db                          | Versio | db                          | Version        |
  | --------------------------- | ------ | --------------------------- | -------------- |
  | NCBI RefSeq viral sequences | latest | ICTV Master Species List    | 2021 v3        |
  | VOG HMM                     | VOG 97 | IMGVR high-quality vOTU rep | v4.1           |
  | iPHoP db                    | latest | GTDB-Tk db                  | release 207 v2 |
  | VirSorter2 db               | latest | DVF db (models)             | v1.0           |

- **Running setting**

  (***note***: here only shows some specific settings; default settings or user-adjustable settings are not included)

  | Command                                          | option                    | explanation                                                  | setting                           |
  | ------------------------------------------------ | ------------------------- | ------------------------------------------------------------ | --------------------------------- |
  | run_virsorter2_1st                               | --keep-original-seq       | keep the original sequences instead of trimmed               | N/A                               |
  |                                                  | --min-score               | minimal score to be identified as viral                      | 0.5                               |
  | run_virsorter2_2nd                               | --seqname-suffix-off      | turn off adding suffix to sequence names                     | N/A                               |
  |                                                  | --viral-gene-enrich-off   | turn off the requirement of more viral than cellular genes for calling full sequence viral | N/A                               |
  |                                                  | --prep-for-dramv          | turn on generating viral seqfile and viral-affi-contigs.tab for DRAMv | N/A                               |
  |                                                  | --min-score               | minimal score to be identified as viral                      | 0.5                               |
  | filter_sorted_bam (coverm)                       | --min-read-aligned-length | exclude reads with smaller numbers of aligned bases          | 50 (short reads) 500 (long reads) |
  | run_vrhyme                                       | --red                     | maximum number of redundant proteins per bin                 | 2                                 |
  | run_vcontact2                                    | --rel-mode                | method to use to create the protein-protein similarity edge file | Diamond                           |
  |                                                  | --db                      | select a reference database to de novo cluster proteins against | None                              |
  |                                                  | --pcs-mode                | whether to use ClusterONE or MCL for Protein Cluster (PC) generation | MCL                               |
  |                                                  | --vcs-mode                | whether to use ClusterONE or MCL for Viral Cluster (VC) generation | ClusterONE                        |
  | run_drep (dRep dereplicate)                      | --ignoreGenomeQuality     | don't run checkM or do any quality filtering                 | N/A                               |
  |                                                  | -pa                       | ANI threshold to form primary (MASH) clusters                | 0.8                               |
  |                                                  | -sa                       | ANI threshold to form secondary clusters                     | 0.95                              |
  |                                                  | -nc                       | minmum level of overlap between genomes when doing secondary comparisons | 0.85                              |
  |                                                  | -comW                     | completeness weight (SCORING CRITERIA)                       | 0                                 |
  |                                                  | -conW                     | contamination weight (SCORING CRITERIA)                      | 0                                 |
  |                                                  | -strW                     | strain heterogeneity weight (SCORING CRITERIA)               | 0                                 |
  |                                                  | -N50W                     | weight of log(genome N50) (SCORING CRITERIA)                 | 0                                 |
  |                                                  | -sizeW                    | weight of log(genome size) (SCORING CRITERIA)                | 1                                 |
  |                                                  | -centW                    | weight of (centrality - S_ani) (SCORING CRITERIA)            | 0                                 |
  | run_diamond_to_RefSeq_viral_protein_db (diamond) | --evalue                  | maximum e-value to report alignments                         | 0.00001                           |
  |                                                  | --query-cover             | minimum query cover% to report an alignment                  | 50                                |
  |                                                  | --subject-cover           | minimum subject cover% to report an alignment                | 50                                |
  |                                                  | -k                        | maximum number of target sequences to report alignments for  | 10000                             |
  | run_hmmsearch_to_marker_VOG_HMM_db (hmmsearch)   | -E                        | report sequences <= this E-value threshold in output         | 0.01                              |
  | run_iphop (iphop predict)                        | --no_qc                   | bypass the automated QC that filters out input sequences with > 10% Ns or with characters other than ATCGN | N/A                               |

  

______
## Running ViWrap <a name="run"></a>

#### ViWrap tasks <a name="flags"></a>

- `run`: Run the full wrapper for identifying, classifying, and characterizing virus genomes from metagenomes

  ***Note***: Before running, please check the header of each input metagenome sequence. Since in the downstream analysis we will add `vRhyme_bin_{i}__` or `vRhyme_unbinned_{i}__` in front of each sequence header, the double understores `__` will be used as the internal splitter. It is suggested to avoid using `__` in your input metagenome, or replace it with other markers.

  ```bash
  # Usage: ViWrap run --input_metagenome <input metagenome assemblies> --input_reads <input metagenomic reads> --out_dir <output directory> [options]
  
  # Example 1 (minimal required inputs):
  ViWrap run --input_metagenome /path/to/Lake_01_assemblies.fasta \
             --input_reads /path/to/Lake_01_T1_1.fastq,/path/to/Lake_01_T1_2.fastq,/path/to/Lake_01_T2_1.fastq,/path/to/Lake_01_T2_2.fastq \
             --out_dir ./ViWrap_Lake_01_outdir \
             --db_dir /path/to/ViWrap_db \
             --identify_method vb \
             --conda_env_dir /path/to/ViWrap_conda_environments
                  
  # Example 2 (with custom MAGs from the same metagenome provided for further host prediction)
  ViWrap run --input_metagenome /path/to/Lake_01_assemblies.fasta \
             --input_reads /path/to/Lake_01_T1_1.fastq,/path/to/Lake_01_T1_2.fastq,/path/to/Lake_01_T2_1.fastq,/path/to/Lake_01_T2_2.fastq \
             --out_dir ./ViWrap_Lake_01_outdir \
             --db_dir /path/to/ViWrap_db \
             --identify_method vs \
             --conda_env_dir /path/to/ViWrap_conda_environments \
             --threads 30 \
             --virome \
             --input_length_limit 5000 \
             --custom_MAGs_dir /path/to/custom_MAGs_dir 
             --iPHoP_db_custom /path/to/iPHoP_db_custom
  ```

- `download`: Download and setup the ViWrap database

  ```bash
  # Note: requires wget, tar, and gzip to be installed
  	
  # Usage:
  ViWrap download --db_dir <output directory for the database>  --conda_env_dir <conda env dir>
  
  # Example:
  ViWrap download --db_dir /path/to/ViWrap_db  --conda_env_dir /path/to/ViWrap_conda_environments
  ```

- `set_up_env`: Set up the conda environments for all scripts 

  ```bash
  # Usage: 
  ViWrap set_up_env --conda_env_dir /path/to/ViWrap_conda_environments  
  ```

- `clean`: Clean redundant information in each result directory

  ```bash
  # Usage:
  ViWrap clean --out_dir ./ViWrap_Lake_01_outdir
      
  # if you have set the "custom_MAGs_dir" option, you can do like this:
  ViWrap clean --out_dir ./ViWrap_Lake_01_outdir --custom_MAGs_dir /path/to/custom_MAGs_dir 
  ```



#### Flag explanations

* `--input_metagenome/-i`: (required) input metagenome assembly. It can be a metagenome or entire virome assembly. The extension of the input nucleotide sequence should be ".fasta".

* `--input_reads/-r`: (required) input metagenomic reads. The input paired reads should be  "forward_1.fastq or forward_R1.fastq" and "reverse_2.fastq or reverse_R2.fastq" connected by ",". Multiple paired reads can be provided at the same time for one metagenome assembly, connected by ",".  For example: `-r /path/to/Lake_01_T1_1.fastq,/path/to/Lake_01_T1_2.fastq,/path/to/Lake_01_T2_1.fastq,/path/to/Lake_01_T2_2.fastq`  Note that the extension of the input reads should be ".fastq" or ".fastq.gz".

* `--input_reads_type/-r`: input metagenomic reads type. The default is "illumina". If you are using long reads, you will need to assign: pacbio - PacBio CLR reads, pacbio_hifi - PacBio HiFi/CCS reads, pacbio_asm20 - PacBio HiFi/CCS reads (asm20), nanopore - Oxford Nanopore reads. We use minimap2 to run long reads

* `--reads_mapping_identity_cutoff/-id`: reads mapping identity cutoff. The default is "0.97". "0.97" is suitable for all illumina reads and also suitable for PacBio Sequel II (HiFi) or Nanopore PromethION Q20+ reads. For other PacBio or Nanopore reads with high error rate, the id cutoff is suggested to be 1 - error rate, i.e., 0.85 (if the error rate is 15%)

* `--out_dir/-o`: (required) output directory to deposit all results (default = ./ViWrap_outdir) output folder to deposit all results. ViWrap will exit if the folder already exists.

* `--db_dir/-d`: (required) database directory (default = $current_dir/ViWrap_db).

* `--identify_method`: (required) the virus identifying method to choose: vb - VIBRANT; vs - VirSorter2 and CheckV; dvf - DeepVirFinder; vb-vs - Use VIBRANT and VirSorter2 to get the overlapped viruses (default); vb-vs-dvf - Use all these three methods and get the overlapped viruses. "vb-vs" is recommended by us since overlapped virus identification will provide more confident results. "vb-vs-dvf" would be too stringent to provide comprehensive virus identification results.

* `--conda_env_dir`: (required) the directory where you put your conda environment files. It is the parent directory that contains all the conda environment folders.

* `--threads/-t`: number of threads (default = 10).

* `--virome/-v`: edit VIBRANT's sensitivity if the input dataset is a virome. It is suggested to use it if you know that the input assembly is virome or metagenome. 

* `--input_length_limit`: length in basepairs to limit input sequences (default=2000, can increase but not decrease); 2000 at least suggested for VIBRANT (vb)-based pipeline, 5000 at least suggested for VirSorter2 (vs)-based pipeline.

* `--custom_MAGs_dir`: custom MAGs dir that contains only *.fasta files for MAGs reconstructed from the same metagenome, this will be used in iPHoP for further host prediction; note that it should be the absolute address path.

* `--iPHoP_db_custom`: custom iPHoP db that will be generated by your input custom MAGs.

* `--iPHoP_db_custom_pre`: custom iPHoP db that has been made from the previous run, this will be used in iPHoP for host prediction by custom db; note that it should be the absolute address path.

  ***Note that***:
  
  (1) if you use `--custom_MAGs_dir`, you will also need to use `--iPHoP_db_custom` to assign the output place of the generated custom iPHoP db.
  
  (2) If you use `--iPHoP_db_custom_pre`, you will not need to use `--custom_MAGs_dir` or `iPHoP_db_custom`. Simply assign the custom iPHoP db that has been made from the previous run.
  
  

#### Test run

```bash
ViWrap run --input_metagenome /path/to/example_data/Guaymas_scaffolds_min1000.subset.fasta \
           --input_reads /path/to/example_data/reads_1.fastq.gz,/path/to/example_data/reads_2.fastq.gz \
           --out_dir ./test_outdir \
           --db_dir /path/to/ViWrap_db \
           --identify_method vb \
           --threads 20 \
           --input_length_limit 5000 \           
           --conda_env_dir /path/to/ViWrap_conda_environments
           
# The total running time is about 2 hrs           
```



______
## Output Explanations <a name="out"></a>

#### **All result folders**
- `00_VIBRANT_input_metagenome_stem_name`: the virus identification result (would be "00_VirSorter_input_metagenome_stem_name", "00_DeepVirFinder_input_metagenome_stem_name", "00_VIBRANT_VirSorter_input_metagenome_stem_name", or "00_VIBRANT_VirSorter_DeepVirFinder_input_metagenome_stem_name")
- `01_Mapping_result_outdir`: the reads mapping result
- `02_vRhyme_outdir`: vRhyme binning result
- `03_vConTACT2_outdir`: vConTACT2 classifying result
- `04_Nlinked_viral_gn_dir`: N-linked viral genome as CheckV inputs
- `05_CheckV_outdir`: CheckV result
- `06_dRep_outdir`: dRep clustering result
- `07_iPHoP_outdir`: iPHoP result for host prediction
- `08_ViWrap_summary_outdir`: Summarized results
- `09_Virus_statistics_visualization`: Visualized statistics of viruses
- `ViWrap_run.log`: running log file containing the issued command and time log

#### **Hierarchy** in `08_ViWrap_summary_outdir`
* '>' : folder
* '-' : file 
```shell
> 08_ViWrap_summary_outdir
    - Genus_cluster_info.txt # Virus genus clusters
    - Species_cluster_info.txt # Virus species clusters
    - Host_prediction_to_genome_m90.csv # Host prediction result at genome level
    - Host_prediction_to_genus_m90.csv # Host prediction result at genus level
    - Sample2read_info.txt # Reads counts and bases
    - Tax_classification_result.txt # Virus taxonomy result
    - Virus_annotation_results.txt # Virus annotation result
    > Virus_genomes_files # Contains all fasta, faa, ffn files for virus genomes
         - vRhyme*.fasta
         - vRhyme*.faa
         - vRhyme*.ffn
    - Virus_normalized_abundance.txt # Normalized virus genome abundance (normalized by 100M reads/sample)
    - Virus_raw_abundance.txt # Raw virus genome abundance
    - Virus_summary_info.txt # Summarized property for all virus genomes
    > AMG_statistics # Contains AMG protein info, statistics, and protein sequences
         - AMG_pro2info.txt
         - Gn2amg_statistics.txt
         - AMG_pros.faa
```

#### **Hierarchy** in `09_Virus_statistics_visualization`

* '>' : folder

* '-' : file 

  ```shell
  > 09_Virus_statistics_visualization
      > Result_visualization_inputs 
          - virus_statistics.txt
          - virus_family_relative_abundance.txt
          - KO_ID_relative_abundance.txt # The relative abundance of KO among all KOs (KO assigned to AMG)
          - KO_metabolism_relative_abundance.txt
      > Result_visualization_outputs    
          - virus_statistics.png # the 1st bar-chart
          - virus_family_relative_abundance.png # the 1st pie-chart
          - KO_ID_relative_abundance.png # the 2nd bar-chart
          - KO_metabolism_relative_abundance.png # the 2nd pie-chart
          - virus_statistics.pdf
          - virus_family_relative_abundance.pdf
          - KO_ID_relative_abundance.pdf
          - KO_metabolism_relative_abundance.pdf
  ```

______
#### Notes

- **The default setting for each software**

  Generally, we all use the default settings for identifying, binning, classifying, and predicting viruses, due to that the default settings are mostly suggested by the inventors and we want to make ViWrap be intuitive and easy to use.

- **Using VIBRANT or VirSorter2 to get prophage information**

  If you want to have prophage information in your final results, we will avoid solely using DVF to identify viruses. You can use "vb", ''vs", "vb-vs" (default), or "vb-vs-dvf" identification methods to get final results that will include prophage information.

- **Using "run_wo_reads" will only provide viral scaffold results**

  - Since vRhyme needs scaffold coverage information for virus binning. Using "run_wo_reads" will not give viral binning results; each viral scaffold will be treated as a viral genome in the downstream analysis after viral binning. In the end, it will generate viral abundance results neither. 
  - It will only generate five result folders: `00_VIBRANT_input_metagenome_stem_name`, `01_vContact2_outdir`, `02_CheckV_outdir`, `03_dRep_outdir`, `04_iPHoP_outdir`, and `05_ViWrap_summary_outdir`
  - This method will facilitate identify viruses from both metagenomes and genomes without the usage of metagenomic/genomic reads.

- **Lytic and lysogenic viruses can bin together** (copied from vRhyme GitHub page)

  * Lytic cycle: productive infection that leads to release of viral particles; strictly lytic viruses do not integrate
  
  * Lysogenic cycle: non-productive infection where viral particles are not released; lysogenic viruses integrate and are dormant until entering the lytic cycle (some exceptions); often encode an integrase
  
  * When viewing binned sequences that have "lysogenic/lytic" labels from another software (e.g., VIBRANT, VirSorter) it may seem concerning to see a lytic virus bin with a lysogenic virus. Although this may be a reason to be skeptical of contamination, here are some explanations of what may be occurring:
    * software tools that label lysogenic/lytic can make mistakes or be misled. For example, if a sequence does not encode an integrase or other lysogenic features then VIBRANT will label it as "lytic". If this is a virus genome fragment and another fragment of the same genome (different sequence) contains an integrase then that other fragment will be labeled as "lysogenic". When binning, those two sequence can be place into the same bin by vRhyme to produce an accurate bin with a lytic and lysogenic member.
    * A bin with one or more lytic members and one lysogenic member should not cause concern.
    * A bin with one or more lytic members and one integrated prophage should be examined.
    * A bin with two or more lysogenic members, each encoding an integrase, is likely contamination.
    * A bin with two or more integrated prophages, regardless of integrases, from multiple parent sequences is likely contamination.
    
    
    
    Based on the above description, here, we provided formulas for automatically assigning the virus bins:
    
    ```bash
    # Formulas:
    # Case 1: 1 lysogenic_scaffold + N lytic_scaffold => lysogenic_virus
    # Case 2: 1 integrated_prophage + N lytic_scaffold => split into scaffolds
    # Case 3: N (N >= 2) of lysogenic_scaffold or integrated_prophage + N lytic_scaffold => split into scaffolds
    # Case 4: N lytic_scaffold => lytic_virus
    ```
    
    
    

______
## Contact <a name="contact"></a>
Please contact Zhichao Zhou (zczhou2017@gmail.com or GitHub Issues) with any questions, concerns or comments.  

Thank you for using ViWrap!  



```
___________________________________________
 __ __  ____  __    __  ____    ____  ____  
|  |  ||    ||  |__|  ||    \  /    ||    \ 
|  |  | |  | |  |  |  ||  D  )|  o  ||  o  )
|  |  | |  | |  |  |  ||    / |     ||   _/ 
|  :  | |  | |  `  '  ||    \ |  _  ||  |   
 \   /  |  |  \      / |  .  \|  |  ||  |   
  \_/  |____|  \_/\_/  |__|\_||__|__||__|   
___________________________________________
```



______
## Copyright  

ViWrap
Copyright (C) 2022 Zhichao Zhou

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
