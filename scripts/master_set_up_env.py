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
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-VIBRANT')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-VIBRANT.yml')} -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-VIBRANT/bin')):
        logger.info(f"{time_current} | ViWrap-VIBRANT conda env has been installed")  
    else:
        logger.info(f"{time_current} | ViWrap-VIBRANT conda env path is not present!")
    
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-geNomad')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-geNomad.yml')} -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-geNomad/bin')):
        logger.info(f"{time_current} | ViWrap-geNomad conda env has been installed")  
    else:
        logger.info(f"{time_current} | ViWrap-geNomad conda env path is not present!")    
        
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-vRhyme.yml')} -y >/dev/null 2>&1")
    #os.system(f"wget https://github.com/AnantharamanLab/vRhyme/raw/master/vRhyme/models/vRhyme_machine_model_ET.sav.gz -q")
    #os.system(f"mv vRhyme_machine_model_ET.sav.gz {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models')}")
    #os.system(f"gzip -d {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models/vRhyme_machine_model_ET.sav.gz')}")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/bin')):
        logger.info(f"{time_current} | ViWrap-vRhyme conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-vRhyme conda env path is not present!")
        
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-vContact2')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-vContact2.yml')} -y >/dev/null 2>&1")
    os.system("wget -c https://paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -q --no-check-certificate")
    
    os.system(f"mv cluster_one-1.0.jar {os.path.join(args['conda_env_dir'], 'ViWrap-vContact2/bin')}")    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-vContact2/bin')):
        logger.info(f"{time_current} | ViWrap-vContact2 conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-vContact2 conda env path is not present!")    
    
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-CheckV')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-CheckV.yml')} -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-CheckV/bin')):
        logger.info(f"{time_current} | ViWrap-CheckV conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-CheckV conda env path is not present!")         
    
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-dRep')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-dRep.yml')} -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-dRep/bin')):
        logger.info(f"{time_current} | ViWrap-dRep conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-dRep conda env path is not present!")     
       
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-Tax')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-Tax.yml')} -y >/dev/null 2>&1") 
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-Tax/bin')):
        logger.info(f"{time_current} | ViWrap-Tax conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-Tax conda env path is not present!") 
    
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-iPHoP')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-iPHoP.yml')} -y >/dev/null 2>&1") 
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-iPHoP/bin')):
        logger.info(f"{time_current} | ViWrap-iPHoP conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-iPHoP conda env path is not present!")     
    
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-GTDBTk')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-GTDBTk.yml')} -y >/dev/null 2>&1")
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-GTDBTk/bin')):
        logger.info(f"{time_current} | ViWrap-GTDBTk conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-GTDBTk conda env path is not present!")       
    
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-vs2')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-vs2.yml')} -y >/dev/null 2>&1")     
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-vs2/bin')):
        logger.info(f"{time_current} | ViWrap-vs2 conda env has been installed")
    else:
        logger.info(f"{time_current} | ViWrap-vs2 conda env path is not present!")          
    
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-Mapping')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-Mapping.yml')} -y >/dev/null 2>&1")   
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-Mapping/bin')):
        logger.info(f"{time_current} | ViWrap-Mapping conda env has been installed")    
    else:
        logger.info(f"{time_current} | ViWrap-Mapping conda env path is not present!")          
    
    os.system(f"mamba env create -p {os.path.join(args['conda_env_dir'], 'ViWrap-DVF')} -f {os.path.join(args['root_dir'],'yamls/ViWrap-DVF.yml')} -y >/dev/null 2>&1")   
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    if os.path.exists(os.path.join(args['conda_env_dir'], 'ViWrap-DVF/bin')):
        logger.info(f"{time_current} | ViWrap-DVF conda env has been installed")     
    else:
        logger.info(f"{time_current} | ViWrap-DVF conda env path is not present!")    
    
    end_time = datetime.now().replace(microsecond=0)
    duration = end_time - start_time
    logger.info(f"The total running time is {duration} (in \"hr:min:sec\" format)")      