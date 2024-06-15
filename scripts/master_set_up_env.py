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
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-VIBRANT')} python=3.7 vibrant=1.2.1 scikit-learn=0.21.3 biopython=1.79 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-VIBRANT/bin')):
        logger.info(f"{time_current} | ViWrap-VIBRANT conda env has been installed")  
    else:
        logger.info(f"{time_current} | ViWrap-VIBRANT conda env path is not present!")
    
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-geNomad')} genomad=1.7.4 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-geNomad/bin')):
        logger.info(f"{time_current} | ViWrap-geNomad conda env has been installed")  
    else:
        logger.info(f"{time_current} | ViWrap-geNomad conda env path is not present!")    
        
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme')} python=3.8 networkx=3.1 pandas=1.0.0 numpy=1.20.3 numba=0.56.4 scikit-learn=0.23.0 pysam=0.21.0 samtools=1.17 mash=2.3 mummer=3.23 mmseqs2=14.7e284 prodigal=2.6.3 bowtie2=2.5.1 bwa=0.7.17 vrhyme=1.1.0 -y >/dev/null 2>&1")
    #os.system(f"wget https://github.com/AnantharamanLab/vRhyme/raw/master/vRhyme/models/vRhyme_machine_model_ET.sav.gz -q")
    #os.system(f"mv vRhyme_machine_model_ET.sav.gz {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models')}")
    #os.system(f"gzip -d {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models/vRhyme_machine_model_ET.sav.gz')}")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/bin')):
        logger.info(f"{time_current} | ViWrap-vRhyme conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-vRhyme conda env path is not present!")
        
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-vContact2')} python=3.7 vcontact2=0.11.0 pytables=3.7.0 biopython=1.79 networkx=2.5.1 numpy=1.19.0 pandas=0.25.3 scipy=1.6.1 scikit-learn=0.24.1 psutil=5.9.3 pyparsing=3.0.9 hdf5=1.12.1 clusterone=1.0 mcl=14.137 blast=2.13.0 diamond=2.0.15 -y >/dev/null 2>&1")
    os.system("wget -c https://paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -q --no-check-certificate")
    
    os.system(f"mv cluster_one-1.0.jar {os.path.join(args['conda_env_dir'], 'ViWrap-vContact2/bin')}")    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-vContact2/bin')):
        logger.info(f"{time_current} | ViWrap-vContact2 conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-vContact2 conda env path is not present!")    
    
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-CheckV')} python=3.8 checkv=1.0.1 diamond=2.0.15 hmmer=3.3.2 prodigal=2.6.3 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-CheckV/bin')):
        logger.info(f"{time_current} | ViWrap-CheckV conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-CheckV conda env path is not present!")         
    
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-dRep')} python=3.11.0 drep=3.4.0 mash=1.1 mummer=3.23 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-dRep/bin')):
        logger.info(f"{time_current} | ViWrap-dRep conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-dRep conda env path is not present!")     
       
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-Tax')} python=3.11.0 diamond=2.0.15 hmmer=3.3.2 -y >/dev/null 2>&1") 
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-Tax/bin')):
        logger.info(f"{time_current} | ViWrap-Tax conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-Tax conda env path is not present!") 
    
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-iPHoP')} python=3.8 iphop=1.3.3 perl-data-dumper perl=5.22.0 -y >/dev/null 2>&1") 
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-iPHoP/bin')):
        logger.info(f"{time_current} | ViWrap-iPHoP conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-iPHoP conda env path is not present!")     
    
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-GTDBTk')} gtdbtk=2.3.2 numpy=1.20.0 -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-GTDBTk/bin')):
        logger.info(f"{time_current} | ViWrap-GTDBTk conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-GTDBTk conda env path is not present!")       
    
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-vs2')} python=3.8 virsorter=2.2.4 numpy=1.20.0 -y >/dev/null 2>&1")     
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-vs2/bin')):
        logger.info(f"{time_current} | ViWrap-vs2 conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-vs2 conda env path is not present!")          
    
    os.system(f"mamba create -c bioconda -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-Mapping')} python=3.8 pysam=0.20.0 bowtie2=2.4.5 coverm=0.6.1 pandas=1.5.3 pyfastx=0.8.4 minimap2=2.24 consent=2.2.2 -y >/dev/null 2>&1")   
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-Mapping/bin')):
        logger.info(f"{time_current} | ViWrap-Mapping conda env has been installed")    
    else:
        logger.info(f"{time_current} | ViWrap-Mapping conda env path is not present!")          
    
    os.system(f"mamba create -c hcc -c conda-forge -p {os.path.join(args['conda_env_dir'], 'ViWrap-DVF')} -c hcc deepvirfinder=2020.11.21 -y >/dev/null 2>&1")   
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-DVF/bin')):
        logger.info(f"{time_current} | ViWrap-DVF conda env has been installed")     
    else:
        logger.info(f"{time_current} | ViWrap-DVF conda env path is not present!")    
    
    end_time = datetime.now().replace(microsecond=0)
    duration = end_time - start_time
    logger.info(f"The total running time is {duration} (in \"hr:min:sec\" format)")      