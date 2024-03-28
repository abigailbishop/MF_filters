import os, sys
import numpy as np
from subprocess import call

def ped_collector(
    Data, Station, Run, analyze_blind_dat = False, include_qual_cut = True
):

    if analyze_blind_dat == False:
        print('chunk_ped is not meant for burn sample! try with 100% data!')
        sys.exit(1)

    print('Collecting Ped starts!')

    from tools.ara_utility import size_checker

    blind_type = ''
    if analyze_blind_dat:
        blind_type = '_full'
    ped_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/ped{blind_type}/'
    if not os.path.exists(ped_path):
        os.makedirs(ped_path)

    if include_qual_cut: 
        qual_path = f'{ped_path}ped{blind_type}_qualities_A{Station}_R{Run}.dat'
        if not os.path.exists(qual_path):
            print(f'{qual_path} is not there!!')
            sys.exit(1)

    out_path = f'{ped_path}ped{blind_type}_values_A{Station}_R{Run}.dat'
    del blind_type, ped_path
    
    repeder_dir = os.environ.get('ARA_UTIL_INSTALL_DIR')+'/bin/'
    print(f'repeder dir: {repeder_dir}')
    os.chdir(repeder_dir)
    del repeder_dir

    if include_qual_cut: 
        print("Running repeder with quality cuts:", qual_path)
        repeder_cmd = f"./repeder -d -m 0 -M 4096 -q {qual_path} {Data} {out_path}"
    else: 
        print("Running repeder without the quality cuts file.")
        repeder_cmd = f"./repeder -d -m 0 -M 4096 {Data} {out_path}"
    print(f'excuted cmd: {repeder_cmd}') 
    call(repeder_cmd.split())
    del repeder_cmd

    print('Ped collecting is done!')
    print(f'output is {out_path}.', size_checker(out_path))

    return


