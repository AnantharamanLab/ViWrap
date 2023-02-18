import sys
import os
import argparse
import logging
import scripts
from scripts import module
from scripts import downloadDB
from datetime import datetime
from pathlib import Path
from glob import glob


def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="download")
    parser.add_argument('--db_dir','-d', dest='db_dir', required=False, default=db_path_default, help=f'database directory; default = {db_path_default}')
    parser.add_argument('--conda_env_dir', dest='conda_env_dir', required=True, default='none', help=r'(required) the directory where you put your conda environment files. It is the parent directory that contains all the conda environment folders')
    parser.add_argument('--threads','-t', dest='threads', required=False, default=10, help='number of threads (default = 10)')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)


def set_defaults(args):
    ## Store databases
    args['CheckV_db'] = os.path.join(args['db_dir'],'CheckV_db')
    args['DRAM_db'] = os.path.join(args['db_dir'],'DRAM_db')
    args['GTDB_db'] = os.path.join(args['db_dir'],'GTDB_db')
    args['iPHoP_db'] = os.path.join(args['db_dir'],'iPHoP_db/iPHoP_db')
    args['iPHoP_db_custom'] = os.path.join(args['db_dir'],'iPHoP_db_custom') 
    args['Kofam_db'] = os.path.join(args['db_dir'],'Kofam_db')
    args['Tax_classification_db'] = os.path.join(args['db_dir'],'Tax_classification_db')
    args['VIBRANT_db'] = os.path.join(args['db_dir'],'VIBRANT_db')
    args['VirSorter2_db'] = os.path.join(args['db_dir'],'VirSorter2_db')
    args['DVF_db'] = os.path.join(args['db_dir'],'DVF_db')    


def main(args):
    # Welcome and logger
    print("### Welcome to ViWrap ###\n") 

	## Set up the logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )    
    logger = logging.getLogger(__name__) 

    # Step 1 Pre-check inputs
    start_time = datetime.now().replace(microsecond=0)

    if not os.path.exists(args['conda_env_dir']):
        sys.exit(f"Could not find conda env dirs within {args['conda_env_dir']}") 
    
    if os.path.exists(args['db_dir']):
        sys.exit(f"The db dir of {args['db_dir']} has also ready been set up")
    else:
        os.mkdir(args['db_dir'])
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | Looks like the input conda software is correct")  
    
    set_defaults(args)
    
    # Step 2  Make VIBRANT db
    vibrant_db_dir_absolute_path = os.path.abspath(args['VIBRANT_db'])
    os.system(f"conda run -p {os.path.join(args['conda_env_dir'], 'ViWrap-VIBRANT')} bash {os.path.join(args['conda_env_dir'], 'ViWrap-VIBRANT/bin/download-db.sh')} {vibrant_db_dir_absolute_path}")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | VIBRANT db has been set up")  
    
    
    # Step 3  Make Tax classification db
    os.mkdir(args['Tax_classification_db'])

    ###############################################
    # Part I Download NCBI RefSeq viral protein db#
    ###############################################

    ## Step 3.1 Download NCBI RefSeq viral protein and protein gpff
    scripts.downloadDB.dl_refseq_viral_protein(args['Tax_classification_db'])
    scripts.downloadDB.dl_refseq_viral_protein_gpff(args['Tax_classification_db'])

    ## Step 3.2 Parse to get protein to NCBI taxonomy info
    scripts.downloadDB.parse_gpff(args['Tax_classification_db'])

    ## Step 3.3 Grep NCBI RefSeq viral proteins with taxonomy info
    scripts.downloadDB.grep_NCBI_RefSeq_viral_proteins_w_tax(args['Tax_classification_db'])

    ## Step 3.4 Reformat NCBI tax to ICTV 8-rank tax
    ictv_tax_info = os.path.join(args['root_dir'], 'database/ICTV_Master_Species_List.txt')
    pro2ictv_8_rank_tax = os.path.join(args['Tax_classification_db'], 'pro2ictv_8_rank_tax.txt')
    scripts.downloadDB.reformat_NCBI_tax_to_ICTV_8_rank_tax(args['Tax_classification_db'], ictv_tax_info, pro2ictv_8_rank_tax)

    ## Step 3.5 Make diamond blastp db
    scripts.downloadDB.make_diamond_db(args['Tax_classification_db'])

    ## Step 3.6 Remove useless files
    scripts.downloadDB.remove(args['Tax_classification_db'])
    
    ##########################
    # Part II Download VOG db#
    ##########################
    ## Step 3.7 Parse to get VOG marker list
    vog_marker_table = os.path.join(args['root_dir'], 'database/VOG_marker_table.txt')
    os.system(f"cp {os.path.join(args['root_dir'], 'database/VOG_marker_table.txt')} {os.path.join(args['Tax_classification_db'], 'VOG_marker_table.txt')}")
    vog_marker_list = scripts.downloadDB.get_vog_marker_table(vog_marker_table)

    ## Step 3.8 Download the latest VOG db and pick VOG markers
    scripts.downloadDB.get_marker_vog_hmm(vog_marker_list, args['Tax_classification_db'])
    
    #############################
    # Part III Download IMGVR db#
    #############################
    ## Step 3.9 cp and degzip IMGVR db
    os.system(f"cat {os.path.join(args['root_dir'], 'database/IMGVR_high-quality_phage_vOTU_representatives.tar.gz*')} > {os.path.join(args['Tax_classification_db'], 'IMGVR_high-quality_phage_vOTU_representatives.tar.gz')}")
    os.system(f"tar xzf {os.path.join(args['Tax_classification_db'], 'IMGVR_high-quality_phage_vOTU_representatives.tar.gz')} --directory {args['Tax_classification_db']}")
    os.system(f"mv {os.path.join(args['Tax_classification_db'], 'IMGVR_high-quality_phage_vOTU_representatives/*')} {args['Tax_classification_db']}")
    os.system(f"rm {os.path.join(args['Tax_classification_db'], 'IMGVR_high-quality_phage_vOTU_representatives.tar.gz')}")
    os.system(f"rm -r {os.path.join(args['Tax_classification_db'], 'IMGVR_high-quality_phage_vOTU_representatives')}")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | Tax classification db has been set up")    
    

    # Step 4 Make CheckV db
    os.system(f"conda run -p {os.path.join(args['conda_env_dir'], 'ViWrap-CheckV')} checkv download_database {args['db_dir']} >/dev/null 2>&1")
    os.system(f"mv {os.path.join(args['db_dir'], 'checkv-db-v*')} {args['CheckV_db']}")
    os.system(f"conda run -p {os.path.join(args['conda_env_dir'], 'ViWrap-CheckV')} export CHECKVDB={args['CheckV_db']}")

    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | CheckV db has been set up")  
    
    
    # Step 5 Make iPHoP db
    os.system(f"wget https://portal.nersc.gov/cfs/m342/iphop/db/iPHoP.latest.tar.gz -O {os.path.join(args['db_dir'], 'iPHoP.latest.tar.gz')}") 
    os.mkdir(os.path.join(args['db_dir'], 'iPHoP_db'))
    os.system(f"tar xzf {os.path.join(args['db_dir'], 'iPHoP.latest.tar.gz')} --directory {os.path.join(args['db_dir'], 'iPHoP_db')}")
    os.system(f"mv {os.path.join(args['db_dir'], 'iPHoP_db/*_pub')} {args['iPHoP_db']}")
    os.system(f"rm {os.path.join(args['db_dir'], 'iPHoP.latest.tar.gz')}")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | iPHoP db has been set up")     
    

    # Step 6 Make GTDB-Tk db (v2.1.1) and release 207 v2
    os.system(f"wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz -O {os.path.join(args['db_dir'], 'gtdbtk_r207_v2_data.tar.gz')} --no-check-certificate")
    os.mkdir(args['GTDB_db'])     
    os.system(f"tar xzf {os.path.join(args['db_dir'], 'gtdbtk_r207_v2_data.tar.gz')} --directory {args['GTDB_db']}")
    os.system(f"mv {os.path.join(args['GTDB_db'], 'release207_v2')} {os.path.join(args['GTDB_db'], 'GTDB_db')}")  
    os.system(f"rm {os.path.join(args['db_dir'], 'gtdbtk_r207_v2_data.tar.gz')}")
    os.system(f"conda env config vars set GTDBTK_DATA_PATH={os.path.join(args['GTDB_db'], 'GTDB_db')} -p {os.path.join(args['conda_env_dir'], 'ViWrap-GTDBTk')}")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | GTDB-Tk db has been set up") 
    

    # Step 7 Download VirSorter2 db
    os.system(f"conda run -p {os.path.join(args['conda_env_dir'], 'ViWrap-vs2')} virsorter setup -d {args['VirSorter2_db']} -j {args['threads']} >/dev/null 2>&1")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | VirSorter2 db has been set up")     
    
    # Step 8 Download DVF db
    os.system(f"git clone https://github.com/jessieren/DeepVirFinder.git {os.path.join(args['db_dir'], 'DVF_db_tmp')}")
    os.system(f"mv {os.path.join(args['db_dir'], 'DVF_db_tmp/models')} {args['DVF_db']}")
    os.system(f"rm -rf {os.path.join(args['db_dir'], 'DVF_db_tmp')}")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | DVF db has been set up")     
    

    end_time = datetime.now().replace(microsecond=0)
    duration = end_time - start_time
    logger.info(f"The total running time is {duration} (in \"hr:min:sec\" format)")  
        