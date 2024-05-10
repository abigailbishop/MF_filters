def make_shell_script(
    file_name, keys, station, runs, 
    blind_dat=1, condor_run=0, not_override=0, 
    l2_data=0, no_tqdm=0, include_qual_cut=True
): 

    shebang = "#!/bin/bash"

    build_environment = [
        "# run the reconstruction script",
        "export HDF5_USE_FILE_LOCKING='FALSE'",
        "source /cvmfs/ara.opensciencegrid.org/trunk/centos7/setup.sh",
        "source /home/abishop/ara/a23/MF_filters/setup.sh",
        "cd /home/abishop/ara/a23/MF_filters/scripts/",
    ]

    debugging_statements = [
        r'printf "\n"',
        'pwd',
        r'printf "\n"',
    ]
    
    commands = []
    squiggles  = "~"*60
    for i, run in enumerate(runs): 
        commands.append(rf'printf "\n\n{squiggles}\n"')
        commands.append(
            fr'printf "\tRunning script for run {run} ({i+1}/{len(runs)})"'
        )
        for key in keys: 
            commands.append(rf'printf "\n{squiggles}\n"')
            commands.append(rf'printf "Running script for {key}\n"')
            qual_type = int(key[-3]) if "qual_cut" in key else 1
            commands.append(
                f"python3 /home/abishop/ara/a23/MF_filters/scripts/script_executor.py "
                f"-k {key} -s {station} -r {run} -b {blind_dat} "
                f"-c {condor_run} -n {not_override} -q {qual_type} "
                f"-t {no_tqdm} -l {l2_data} -qc {include_qual_cut} "
            )
            commands.append(rf'printf "Done!\n\n\n"')
            commands.append(rf'printf "\n"')
        commands.append(rf'printf "\n"')
        commands.append("")

    file = open(file_name, "w")
    file.write(shebang+"\n\n")
    for block in [build_environment, debugging_statements, commands]:
        for line in block: 
            file.write(line+"\n") 
        file.write("\n\n")
    file.close()

    return None

def main():

    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runlist', type=str,
        help="Absolute path to runlist")
    parser.add_argument('output_dir', type=str,
        help="Absolute path to the directory where new shell scripts will be saved.")
    parser.add_argument('station', type=int,
        help="Number of the ARA station")
    parser.add_argument('-k', '--keys', type=str, nargs="*",
        help="key corresponding to the script called by script_executor.py")
    parser.add_argument('--first_run', type=int, default = 0,
        help="First run being analyzed")
    parser.add_argument('--last_run', type=int, default = 100000,
        help="Last run being analyzed")
    parser.add_argument('--runs_per_job', type=int, default = 10,
        help="Number of runs to execute the script for per job.")
    parser.add_argument('--blind_dat', type=int, default=1)
    parser.add_argument('--condor_run', type=int, default=0)
    parser.add_argument('--not_override', type=int, default=0)
    parser.add_argument('--l2_data', type=int, default=0)
    parser.add_argument('--no_tqdm', type=int, default=0)
    parser.add_argument('--include_qual_cut', type=bool, default=True)
    parser.add_argument('--missing_from_dir', default="", type=str)
    parser.add_argument('--full_pipeline', action="store_true")

    args = parser.parse_args()

    if args.missing_from_dir: 
        from find_missing_outputs import find_missing_outputs
        input_runs = find_missing_outputs(
            args.runlist, args.missing_from_dir, 
            first_run=args.first_run, last_run=args.last_run)
    else: 
        from find_missing_outputs import get_runs_from_runlist
        input_runs = get_runs_from_runlist(
            args.runlist, args.first_run, args.last_run)
    
    split_runs = [
        input_runs[ i*args.runs_per_job : (i+1)*args.runs_per_job ]
        if (i*args.runs_per_job) < len(input_runs)
        else input_runs[ i*args.runs_per_job : ]
        for i in range( (len(input_runs)//args.runs_per_job) + 1  )
    ]
    print(f"Making {len(split_runs)} shells scripts for {len(input_runs)} jobs.")

    if args.full_pipeline: 
        keys = [
            'sub_info', 'qual_cut_1st', 'ped', 'baseline', 
            'cw_flag', 'cw_band', 'cw_ratio', 'qual_cut_2nd',
            'ped', 'rayl', 'rayl_nb', 'snr', 
            'reco_ele_lite', 'rpr', 'vertex', 'mf', 'qual_cut_3rd',
            'sub_info_burn'
        ]
    else: 
        keys = args.keys

    for i, runs in enumerate(split_runs): 
        make_shell_script(
            f"{args.output_dir}/job{i}.sh", keys, args.station, runs,
            blind_dat=args.blind_dat, condor_run=args.condor_run, 
            not_override=args.not_override, 
            l2_data=args.l2_data, no_tqdm=args.no_tqdm, 
            include_qual_cut=args.include_qual_cut
        )


if __name__=="__main__":
    main()