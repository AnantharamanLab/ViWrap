import sys
import os
import argparse
import logging
from datetime import datetime


def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="set_up_env")
    parser.add_argument('--conda', dest='conda', required=True, default='none', help='conda software: miniconda3 or anaconda3')
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

    if not args['conda']:
        sys.exit(f"Could not find the input of conda software, for example: miniconda3 or anaconda3")
    elif not os.path.exists(os.path.join(os.path.expanduser('~'), args['conda'], 'envs')):
        sys.exit(f"Could not find the conda envs folder of {args['conda']}")

    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | Looks like the input conda software is correct")
             
    # Step 2 Install conda env
    os.system("mamba create -c bioconda -n ViWrap-VIBRANT python=3.7 vibrant==1.2.1 scikit-learn==0.21.3 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-VIBRANT conda env has been installed")
    
    os.system("mamba create -c bioconda -n ViWrap-vRhyme python=3.8 networkx pandas numpy numba scikit-learn pysam samtools mash mummer mmseqs2 prodigal bowtie2 bwa vrhyme -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-vRhyme conda env has been installed")
    
    os.system("mamba create -n ViWrap-vContact2 python=3.8 vcontact2==0.11.3 scipy=1.6.0 mcl blast diamond clusterone -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-vContact2 conda env has been installed")
    
    os.system("mamba create -c bioconda -n ViWrap-CheckV python=3 checkv diamond hmmer prodigal -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-CheckV conda env has been installed")
    
    os.system("mamba create -c bioconda -n ViWrap-dRep python=3 drep mash mummer -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-dRep conda env has been installed")
    
    os.system("mamba create -c bioconda -n ViWrap-Tax python=3 diamond==2.0.15 hmmer -y >/dev/null 2>&1") 
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-Tax conda env has been installed")
    
    os.system("mamba create -c bioconda -n ViWrap-iPHoP python=3.8 iphop==1.1.0 -y >/dev/null 2>&1") 
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-iPHoP conda env has been installed")
    
    os.system("mamba create -c bioconda -n ViWrap-GTDBTk gtdbtk==1.6.0 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-GTDBTk conda env has been installed")
    
    os.system("mamba create -c conda-forge -c bioconda -n ViWrap-vs2 virsorter==2 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-vs2 conda env has been installed")
    
    os.system("mamba create -c bioconda -n ViWrap-DRAM dram==1.3.5 -y >/dev/null 2>&1")   
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | ViWrap-DRAM conda env has been installed")
    
    # Step 3 Move files
    os.system("wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -q")
    os.system(f"mv cluster_one-1.0.jar {os.path.join(os.path.expanduser('~'), args['conda'], 'envs/ViWrap-vContact2/bin')}")

    os.system(f"wget https://github.com/AnantharamanLab/vRhyme/raw/master/vRhyme/models/vRhyme_machine_model_ET.sav.gz -q")
    os.system(f"mv vRhyme_machine_model_ET.sav.gz {os.path.join(os.path.expanduser('~'), args['conda'], 'envs/ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models')}")
    os.system(f"gzip -d {os.path.join(os.path.expanduser('~'), args['conda'], 'envs/ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models/vRhyme_machine_model_ET.sav.gz')}")
    
    # Step 4 Check installed conda env
    conda_env_list_call_out = os.popen('conda env list').read()
    conda_env_list_call_out_list = conda_env_list_call_out.split('\n')
    ViWrap_env_map = {} # env => addr
    for line in conda_env_list_call_out_list:
        if line.startswith('ViWrap-'):
            tmp = line.split()
            env, addr = tmp[0], tmp[1]
            ViWrap_env_map[env] = addr
            
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

    if os.path.exists(ViWrap_env_map['ViWrap-DRAM']):
        logger.info("ViWrap-DRAM conda env path has been checked")
    else:
        logger.info("ViWrap-DRAM conda env path is not present!")  

    end_time = datetime.now().replace(microsecond=0)
    duration = end_time - start_time
    logger.info(f"The total running time is {duration} (in \"hr:min:sec\" format)")      
    
    
    
    

