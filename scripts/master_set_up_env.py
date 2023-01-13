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

def install(prefix, env, channel, pkgs, logger):
    env_path = os.path.join(prefix, env)
    # will silently fail since redirecting stdout/stderr to /dev/null
    cmd = f"mamba create -c conda-forge -c {channel} -p {env_path} {pkgs} -y > /dev/null 2>&1"
    os.system(cmd)
    
    # thus need to check install IMMEDIATELY
    time = f"[{str(datetime.now().replace(microsecond=0))}]"
    if not os.path.exists(os.path.join(env_path, "bin")):
        logger.info(f"{time} | {env} conda env was not installed correctly")
        raise RuntimeError(f"{env} could not be installed correctly using this conda command:\n{cmd}\n")
    else:
        logger.info(f"{time} | {env} conda env path has been installed")

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
    env_pkgs = {
        "ViWrap-VIBRANT": (
            "python=3.7 vibrant==1.2.1 scikit-learn==0.21.3 biopython", 
            "bioconda"
        ),
        "ViWrap-vRhyme": (
            "python=3.8 networkx pandas==1.0.0 numpy==1.20 numba scikit-learn==0.23.0 pysam samtools mash mummer mmseqs2 prodigal bowtie2 bwa vrhyme=1.1.0", 
            "bioconda"
        ),
        "ViWrap-vContact2": (
            "python=3.7 vcontact2=0.11.0 pytables biopython networkx numpy=1.19.0 pandas=0.25.3 scipy==1.6.1 scikit-learn==0.24.1 psutil pyparsing hdf5 clusterone mcl blast diamond", 
            "bioconda"
        ),
        "ViWrap-CheckV": (
            "python=3.8 checkv diamond hmmer prodigal", 
            "bioconda"),
        "ViWrap-dRep": (
            "python=3 drep mash mummer", 
            "bioconda"
        ),
        "ViWrap-Tax": (
            "python=3 diamond==2.0.15 hmmer", 
            "bioconda"
        ),
        "ViWrap-iPHoP": (
            "python=3.8 iphop==1.1.0", 
            "bioconda"
        ),
        "ViWrap-GTDBTk": (
            "gtdbtk==2.1.1 numpy=1.20.0", 
            "bioconda"
        ),
        "ViWrap-vs2": (
            "python=3.8 virsorter==2.2.3 numpy=1.20.0", 
            "bioconda"
        ),
        "ViWrap-Mapping": (
            "python=3.8 pysam bowtie2==2.4.5 coverm pandas pyfastx minimap2 consent", 
            "bioconda"
        ),
        "ViWrap-DVF": (
            "deepvirfinder", 
            "hcc"
        ),
    }

    for env, (pkgs, channel) in env_pkgs.items():
        install(args['conda_env_dir'], env, channel, pkgs, logger)

    # additional step for virsorter2
    # Add "  - numpy=1.20.0" to the conda env yaml file
    vs2_yaml_file = os.path.join(args['conda_env_dir'], 'ViWrap-vs2', 'lib/python3.8/site-packages/virsorter/envs/vs2.yaml')
    with open(vs2_yaml_file, 'a') as f:
        f.write('  - numpy=1.20.0\n')    
    
    
    # Step 3 Move files
    os.system("wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -q")
    os.system(f"mv cluster_one-1.0.jar {os.path.join(args['conda_env_dir'], 'ViWrap-vContact2/bin')}")

    os.system(f"wget https://github.com/AnantharamanLab/vRhyme/raw/master/vRhyme/models/vRhyme_machine_model_ET.sav.gz -q")
    os.system(f"mv vRhyme_machine_model_ET.sav.gz {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models')}")
    os.system(f"gzip -d {os.path.join(args['conda_env_dir'], 'ViWrap-vRhyme/lib/python3.8/site-packages/vRhyme/models/vRhyme_machine_model_ET.sav.gz')}")          

    end_time = datetime.now().replace(microsecond=0)
    duration = end_time - start_time
    logger.info(f"The total running time is {duration} (in \"hr:min:sec\" format)")      
   
    
    
    

