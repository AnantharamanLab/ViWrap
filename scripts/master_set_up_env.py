import sys
import os
import argparse
import logging
from datetime import datetime


def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="set_up_env")
    parser.add_argument('--conda_env_dir', dest='conda_env_dir', required=True, default='none', help=r'(required) the directory where you put your conda environment files. It is the parent directory that contains all the conda environment folders')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)
    
def main(args):
    # Welcome and logger
    print("### Set up conda env ###\n") 
    
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

    if not args['conda_env_dir']:
        sys.exit(f"Could not find the input parameter of conda environments directory")
    elif not os.path.exists(args['conda_env_dir']):
        os.mkdir(args['conda_env_dir'])    

    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | Looks like the input parameter is correct")
             
    # Step 2 Install conda env
    
    os.system(f"mamba create -c bioconda -p {os.path.join(args['conda_env_dir'], 'ViWrap-VIBRANT')} python=3.7 vibrant==1.2.1 scikit-learn==0.21.3 biopython -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-VIBRANT conda env has been installed")
    
    os.system(f"mamba create -c bioconda -p {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme')} python=3.8 networkx pandas numpy numba scikit-learn pysam samtools mash mummer mmseqs2 prodigal bowtie2 bwa vrhyme -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-vRhyme conda env has been installed")
    
    os.system(f"mamba create -p {os.path.join(args['conda_env_dir'], 'ViWrap-vContact2')} python=3.8 vcontact2==0.11.3 scipy=1.6.0 mcl blast diamond clusterone -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-vContact2 conda env has been installed")
    
    os.system(f"mamba create -c bioconda -p {os.path.join(args['conda_env_dir'], 'ViWrap-CheckV')} python=3 checkv diamond hmmer prodigal -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-CheckV conda env has been installed")
    
    os.system(f"mamba create -c bioconda -p {os.path.join(args['conda_env_dir'], 'ViWrap-dRep')} python=3 drep mash mummer -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-dRep conda env has been installed")
    
    os.system(f"mamba create -c bioconda -p {os.path.join(args['conda_env_dir'], 'ViWrap-Tax')} python=3 diamond==2.0.15 hmmer -y >/dev/null 2>&1") 
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-Tax conda env has been installed")
    
    os.system(f"mamba create -c bioconda -p {os.path.join(args['conda_env_dir'], 'ViWrap-iPHoP')} python=3.8 iphop==1.1.0 -y >/dev/null 2>&1") 
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-iPHoP conda env has been installed")
    
    os.system(f"mamba create -c bioconda -p {os.path.join(args['conda_env_dir'], 'ViWrap-GTDBTk')} gtdbtk==1.6.0 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-GTDBTk conda env has been installed")
    
    os.system(f"mamba create -c conda-forge -c bioconda -p {os.path.join(args['conda_env_dir'], 'ViWrap-vs2')} virsorter==2.2.3 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-vs2 conda env has been installed")
    
    os.system(f"mamba create -c bioconda -p {os.path.join(args['conda_env_dir'], 'ViWrap-Mapping')} python=3 pysam bowtie2==2.4.5 coverm pandas pyfastx -y >/dev/null 2>&1")   
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-Mapping conda env has been installed")    
    
    os.system(f"mamba create -c hcc -p {os.path.join(args['conda_env_dir'], 'ViWrap-DVF')} python=3 deepvirfinder -y >/dev/null 2>&1")   
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-DVF conda env has been installed")     
    
    
    # Step 3 Move files
    os.system("wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -q")
    os.system(f"mv cluster_one-1.0.jar {os.path.join(args['conda_env_dir'], 'ViWrap-vContact2/bin')}")

    os.system(f"wget https://github.com/AnantharamanLab/vRhyme/raw/master/vRhyme/models/vRhyme_machine_model_ET.sav.gz -q")
    os.system(f"mv vRhyme_machine_model_ET.sav.gz {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models')}")
    os.system(f"gzip -d {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models/vRhyme_machine_model_ET.sav.gz')}")
    
    # Step 4 Check installed conda env
    conda_env_list_call_out = os.popen('conda env list').read()
    conda_env_list_call_out_list = conda_env_list_call_out.split('\n')
    ViWrap_env_map = {} # env to addr 
    for line in conda_env_list_call_out_list:
        if 'ViWrap-' in line:
            line = line.rstrip('\n')
            line = line.strip()
            env = line.rsplit('/', 1)[1]
            ViWrap_env_map[env] = line
            
    if os.path.exists(ViWrap_env_map['ViWrap-VIBRANT']):
        logger.info("ViWrap-VIBRANT conda env path has been checked")
    else:
        logger.info("ViWrap-VIBRANT conda env path is not present!")
        
    if os.path.exists(ViWrap_env_map['ViWrap-vRhyme']):
        logger.info("ViWrap-vRhyme conda env path has been checked")
    else:
        logger.info("ViWrap-vRhyme conda env path is not present!")

    if os.path.exists(ViWrap_env_map['ViWrap-vContact2']):
        logger.info("ViWrap-vContact2 conda env path has been checked")
    else:
        logger.info("ViWrap-vContact2 conda env path is not present!")

    if os.path.exists(ViWrap_env_map['ViWrap-CheckV']):
        logger.info("ViWrap-CheckV conda env path has been checked")
    else:
        logger.info("ViWrap-CheckV conda env path is not present!")

    if os.path.exists(ViWrap_env_map['ViWrap-dRep']):
        logger.info("ViWrap-dRep conda env path has been checked")
    else:
        logger.info("ViWrap-dRep conda env path is not present!")

    if os.path.exists(ViWrap_env_map['ViWrap-Tax']):
        logger.info("ViWrap-Tax conda env path has been checked")
    else:
        logger.info("ViWrap-Tax conda env path is not present!")

    if os.path.exists(ViWrap_env_map['ViWrap-iPHoP']):
        logger.info("ViWrap-iPHoP conda env path has been checked")
    else:
        logger.info("ViWrap-iPHoP conda env path is not present!")

    if os.path.exists(ViWrap_env_map['ViWrap-GTDBTk']):
        logger.info("ViWrap-GTDBTk conda env path has been checked")
    else:
        logger.info("ViWrap-GTDBTk conda env path is not present!")

    if os.path.exists(ViWrap_env_map['ViWrap-vs2']):
        logger.info("ViWrap-vs2 conda env path has been checked")
    else:
        logger.info("ViWrap-vs2 conda env path is not present!")  

    if os.path.exists(ViWrap_env_map['ViWrap-DVF']):
        logger.info("ViWrap-DVF conda env path has been checked")
    else:
        logger.info("ViWrap-DVF conda env path is not present!")  
        
    if os.path.exists(ViWrap_env_map['ViWrap-Mapping']):
        logger.info("ViWrap-Mapping conda env path has been checked")
    else:
        logger.info("ViWrap-Mapping conda env path is not present!")          

    end_time = datetime.now().replace(microsecond=0)
    duration = end_time - start_time
    logger.info(f"The total running time is {duration} (in \"hr:min:sec\" format)")      
    
    
    
    

