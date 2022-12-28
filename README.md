<p align="center">
  <img width="350"  src="https://github.com/AnantharamanLab/ViWrap/blob/main/images/ViWrap2.jpg">
</p>

# **ViWrap**
### A wrapper to identify, bin, classify, and predict host-viral relationship for viruses from Metagenomes

```
Oct 2022   
Zhichao Zhou  
zzhou388@wisc.edu & zczhou2017@gmail.com  
Anantharaman Lab
Department of Bacteriology
University of Wisconsin-Madison  
```

## Current Version
ViWrap v1.2.0 

<img src="https://github.com/AnantharamanLab/ViWrap/blob/main/images/ViWrap-Flowchart.jpg">

## Citation
If you find ViWrap useful please consider citing our manuscript on XXX (still in draft):  

```tex
Zhou, Z et al. ViWrap, a pipeline to identify, bin, classify, and predict host-viral relationship for viruses from Metagenomes. YYY. 2023
```

______
## Table of Contents:
1. [Updates](#updates)
    * v1.1.0
    * v1.0.0
2. [Program Description](#program)
3. [Installation](#install)
5. [Running ViWrap](#run)
    * ViWrap tasks
    * Flag explanations
5. [Output Explanations](#out)
10.  [Contact](#contact)

______
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
* Add phylogenetic tree building step for Caudovirales



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

[VIBRANT](https://github.com/AnantharamanLab/VIBRANT): https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0

[VirSorter2](https://github.com/jiarong/VirSorter2): https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y

[DeepVirFinder](https://github.com/jessieren/DeepVirFinder): https://link.springer.com/article/10.1007/s40484-019-0187-4

[vContact2](https://bitbucket.org/MAVERICLab/vcontact2/wiki/Home): https://www.nature.com/articles/s41587-019-0100-8

[vRhyme](https://github.com/AnantharamanLab/vRhyme): https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac341/6584432

[iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/) ([and software within it](https://bitbucket.org/srouxjgi/iphop/src/main/#markdown-header-citation)): https://www.biorxiv.org/content/10.1101/2022.07.28.501908v1


#### **ViWrap Features**
* Identify and annotate viruses by VIBRANT, VirSorter2, or DeepVirFinder
* Map reads onto all input metagenome assemblies to get scaffold depth
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

Since ViWrap has many dependencies to be installed, it would be much easier to set up a conda environment instead of installing all dependencies in the global environment (make sure you have upfront conda installed on your server, i.e., [miniconda3](https://docs.conda.io/en/latest/miniconda.html) or anaconda3; we only suggest to run in version 3.0+ conda). Since ViWrap will use multiple conda environments with considerably large sizes, we strongly suggest placing them elsewhere (for example, here it is `/path/to/ViWrap_conda_environments`) instead of the home address ( `$HOME` ) as normally done.

1. `conda create -c bioconda -p /path/to/ViWrap_conda_environments/ViWrap python=3.8 biopython mamba numpy pandas pyfastx matplotlib seaborn`
2. `conda activate /path/to/ViWrap_conda_environments/ViWrap`

Note: `/path/to/conda_environments` indicates the directory that you will need to use to store all conda environments for ViWrap

#### GitHub installation

1. `git clone https://github.com/AnantharamanLab/ViWrap`
2. `cd ViWrap`
3. `chmod +x ViWrap scripts/*.py # Make all python scripts to be executable`
4. `PATH=`pwd`:$PATH # Add ViWrap to the PATH, so it can be called elsewhere in a terminal`

#### Set up the other conda environments required by ViWrap

`ViWrap set_up_env --conda_env_dir /path/to/ViWrap_conda_environments`

This will take several minutes depending on your current internet speed. 

**Note**: `/path/to/ViWrap_conda_environments` can be set anywhere on your server to contain ViWrap conda environments, it does not necessarily have to be "miniconda3/envs" or "anaconda3/envs"; on contrary, we suggest avoiding these two addresses, it is better to be set to elsewhere that is not within `$HOME`.

ViWrap will use the "-p" or "--prefix" option to specify where to write the environment files:

(Details in https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#specifying-a-location-for-an-environment)

For example,`conda create --prefix /tmp/test-env python=3.8`
will create the environment named `/tmp/test-env` which resides in `/tmp/` instead of the default `.conda`.

The following 11 conda environments will be set up:

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

#### Set up ViWrap database

```shell
ViWrap download --db_dir /path/to/ViWrap_db  --conda_env_dir /path/to/ViWrap_conda_environments
```

`/path/to/ViWrap_db` is the place you store the ViWrap database. Please make sure there is enough space to store the database (~280G at least). It will take ~3-4 hours to set up well depending on your current internet speed. This is kind of tedious, however, you will only need to do this one time.

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

#### Test ViWrap

1. `ViWrap -h`
2. `ViWrap run -h`

______
## Running ViWrap <a name="run"></a>

#### ViWrap tasks <a name="flags"></a>

- `run`: Run the full wrapper for identifying, classifying, and characterizing virus genomes from metagenomes

  ```python
  # Usage: ViWrap run --input_metagenome <input metagenome assemblies> --input_reads <input metagenomic reads> --out_dir <output directory> [options]
  
  # Example 1 (minimal required inputs):
  ViWrap run --input_metageome /path/to/Lake_01_assemblies.fasta \
             --input_reads /path/to/Lake_01_T1_1.fastq,/path/to/Lake_01_T1_2.fastq,/path/to/Lake_01_T2_1.fastq,/path/to/Lake_01_T2_2.fastq \
             --out_dir ./ViWrap_Lake_01_outdir \
             --db_dir /path/to/ViWrap_db \
             --identify_method vb \
             --conda_env_dir /path/to/ViWrap_conda_environments
                  
  # Example 2 (with custom MAGs from the same metagenome provided for further host prediction)
  ViWrap run --input_metageome /path/to/Lake_01_assemblies.fasta \
             --input_reads /path/to/Lake_01_T1_1.fastq,/path/to/Lake_01_T1_2.fastq,/path/to/Lake_01_T2_1.fastq,/path/to/Lake_01_T2_2.fastq \
             --out_dir ./ViWrap_Lake_01_outdir \
             --db_dir /path/to/ViWrap_db \
             --identify_method vs \
             --conda_env_dir /path/to/ViWrap_conda_environments \
             --threads 30 \
             --virome \
             --input_length_limit 2000 \
             --custom_MAGs_dir /path/to/custom_MAGs_dir 
              
  ```

- `download`: Download and setup the ViWrap database

  ```python
  # Note: requires wget, tar, and gzip to be installed
  	
  # Usage:
  ViWrap download --db_dir <output directory for the database>  --conda_env_dir <conda env dir>
  
  # Example:
  ViWrap download --db_dir /path/to/ViWrap_db  --conda_env_dir /path/to/ViWrap_conda_environments
  ```

- `set_up_env`: Set up the conda environments for all scripts 

  ```python
  # Usage: 
  ViWrap set_up_env --conda_env_dir /path/to/ViWrap_conda_environments  
  ```

- `clean`: Clean redundant information in each result directory

  ```python
  # Usage:
  ViWrap clean --out_dir ./ViWrap_Lake_01_outdir
      
  # if you have set the "custom_MAGs_dir" option, you can do like this:
  ViWrap clean --out_dir ./ViWrap_Lake_01_outdir --custom_MAGs_dir /path/to/custom_MAGs_dir 
  ```



#### Flag explanations

* `--input_metagenome/-i`: (required) input metagenome assembly. It can be a metagenome or entire virome assembly. The extension of the input nucleotide sequence should be ".fasta".
* `--input_reads/-r`: (required) input metagenomic reads. The input paired reads should be  "forward_1.fastq or forward_R1.fastq" and "reverse_2.fastq or reverse_R2.fastq" connected by ",". Multiple paired reads can be provided at the same time for one metagenome assembly, connected by ",".  For example: `-r /path/to/Lake_01_T1_1.fastq,/path/to/Lake_01_T1_2.fastq,/path/to/Lake_01_T2_1.fastq,/path/to/Lake_01_T2_2.fastq`  Note that the extension of the input reads should be ".fastq" or ".fastq.gz".

* `--out_dir/-o`: (required) output directory to deposit all results (default = ./ViWrap_outdir) output folder to deposit all results. ViWrap will exit if the folder already exists.
* `--db_dir/-d`: (required) database directory (default = $current_dir/ViWrap_db).
* `--identify_method`: (required) the virus identifying method to choose: vb - VIBRANT; vs - VirSorter2 and CheckV; dvf - DeepVirFinder; vb-vs - Use VIBRANT and VirSorter2 to get the overlapped viruses (default); vb-vs-dvf - Use all these three methods and get the overlapped viruses. "vb-vs" is recommended by us since overlapped virus identification will provide more confident results. "vb-vs-dvf" would be too stringent to provide comprehensive virus identification results.

* `--conda_env_dir`: (required) the directory where you put your conda environment files. It is the parent directory that contains all the conda environment folders.
* `--threads/-t`: number of threads (default = 10).
* `--virome/-v`: edit VIBRANT's sensitivity if the input dataset is a virome. It is suggested to use it if you know that the input assembly is virome or metagenome. 

* `--input_length_limit`: length in basepairs to limit input sequences (default=2000, can increase but not decrease); 2000 at least suggested for VIBRANT (vb)-based and INHERIT (in)-based pipeline, 5000 at least suggested for VirSorter2 (vs)-based pipeline.
* `--custom_MAGs_dir`: custom MAGs dir that contains only *.fasta files for MAGs reconstructed from the same metagenome, this will be used in iPHoP for further host prediction; note that it should be the absolute address path.

______
## Output Explanations <a name="out"></a>

#### **All result folders**
- `00_VIBRANT_input_metageome_stem_name`: the virus identification result (would be "00_VirSorter_input_metageome_stem_name", "00_DeepVirFinder_input_metageome_stem_name", "00_VIBRANT_VirSorter_input_metageome_stem_name", or "00_VIBRANT_VirSorter_DeepVirFinder_input_metageome_stem_name")
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
    > Virus_normalized_abundance.txt # Normalized virus genome abundance (normalized by 100M reads/sample)
    > Virus_raw_abundance.txt # Raw virus genome abundance
    > Virus_summary_info.txt # Summarized property for all virus genomes
```

#### **Hierarchy** in `09_Virus_statistics_visualization`

* '>' : folder

* '-' : file 

  ```shell
  > 09_Virus_statistics_visualization
      > Result_visualization_inputs 
          - virus_statistics.txt
          - virus_family_relative_abundance.txt
          - KO_ID_relative_abundance.txt
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
  - It will only generate five result folders: `00_VIBRANT_input_metageome_stem_name`, `01_vContact2_outdir`, `02_CheckV_outdir`, `03_dRep_outdir`, `04_iPHoP_outdir`, and `05_ViWrap_summary_outdir`
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
