import sys
import os
import argparse
import logging
from pathlib import Path
from datetime import datetime


def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="clean")
    parser.add_argument('--out_dir','-o', dest='out_dir', required=False, default='./ViWrap_outdir', help=r'(required) output directory to deposit all results (default = ./ViWrap_outdir) output folder to deposit all results. ViWrap will exit if the folder already exists')
    parser.add_argument('--custom_MAGs_dir', dest='custom_MAGs_dir', required=False, default='none', help=r'custom MAGs dir that contains only *.fasta files for MAGs reconstructed from the same metagenome, this will be used in iPHoP for host prediction; note that it should be the absolute address path')	
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)

def set_defaults(args):
    ## Store outdirs 
    args['mapping_outdir'] = os.path.join(args['out_dir'],'01_Mapping_result_outdir')
    args['vrhyme_outdir'] = os.path.join(args['out_dir'],'02_vRhyme_outdir')
    args['vcontact2_outdir'] = os.path.join(args['out_dir'],'03_vContact2_outdir')
    args['nlinked_viral_gn_dir'] = os.path.join(args['out_dir'],'04_Nlinked_viral_gn_dir')
    args['checkv_outdir'] = os.path.join(args['out_dir'],'05_CheckV_outdir')
    args['drep_outdir'] = os.path.join(args['out_dir'],'06_dRep_outdir')
    args['iphop_outdir'] = os.path.join(args['out_dir'],'07_iPHoP_outdir')
    args['iphop_custom_outdir'] = os.path.join(args['out_dir'],'07_iPHoP_outdir/iPHoP_outdir_custom_MAGs')
    args['viwrap_summary_outdir'] = os.path.join(args['out_dir'],'08_ViWrap_summary_outdir')
    
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

    if not os.path.exists(args['out_dir']):
        sys.exit(f"Make sure you give the correct path of ViWrap outdir")       

    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | Looks like the input parameters are correct")  
    
    set_defaults(args)

    # Step 2 Clean 01_Mapping_result_outdir
    os.system(f"rm {os.path.join(args['mapping_outdir'], '*bowtie2_idx*')}")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | 01_Mapping_result_outdir has been cleaned")  
    
    
    # Step 3 Clean 02_vRhyme_outdir    
    os.system(f"rm {os.path.join(args['vrhyme_outdir'], 'combined_pro2viral_gn_map.csv')}")
    os.system(f"rm {os.path.join(args['vrhyme_outdir'], 'combined_viral_faa.faa')}")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | 02_vRhyme_outdir has been cleaned") 
    
    
    # Step 4 Clean 03_vContact2_outdir
    os.system(f"rm {os.path.join(args['vcontact2_outdir'], 'combined_viral_faa.*')}")
    os.system(f"rm {os.path.join(args['vcontact2_outdir'], 'modules.ntwk')}")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | 03_vContact2_outdir has been cleaned")
    

    # Step 5 Clean 07_iPHoP_outdir
    os.system(f"rm -r {os.path.join(args['iphop_outdir'], 'Wdir')}")
    if args['custom_MAGs_dir'] != 'none':
        os.system(f"rm -r {os.path.join(args['iphop_custom_outdir'], 'Wdir')}")
    
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | 07_iPHoP_outdir has been cleaned")

    end_time = datetime.now().replace(microsecond=0)
    duration = end_time - start_time
    logger.info(f"The total running time is {duration} (in \"hr:min:sec\" format)")