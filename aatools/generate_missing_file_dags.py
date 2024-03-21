JOB job_ARA_S2_R19819 ARA_job.sub
VARS job_ARA_S2_R19819 station="2" run="19819"

def generate_dag(runlist, save_dir, station): 
    from numpy import savetxt
    output = []
    for run in runlist: 
        output.append(f"JOB job_ARA_S2_R{run} ARA_job.sub")
        output.append(f'VARS job_ARA_S2_R{run} station="{station}" run="{run}"')
    savetxt(save_dir+f"A{station}_missing.dag", output)


def main():

    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runlist', type=str,
        help="Absolute path to runlist")
    parser.add_argument('output_dir', type=str,
        help="Absolute path to output files.")
    parser.add_argument('save_dir', type=str,
        help="Absolute path to where the new dag file should be saved")
    parser.add_argument('station', type=int,
        help="Number of the ARA station")
    parser.add_argument('--first_run', type=int, default = 0,
        help="First run being analyzed")
    parser.add_argument('--last_run', type=int, default = 100000,
        help="Last run being analyzed")
    args = parser.parse_args()

    from find_missing_outputs import find_missing_outputs
    missing_runs = find_missing_outputs(
        args.runlist, args.output_dir, args.first_run, args.last_run)
    generate_dag(missing_runs, args.savedir, args.station)

    return

if __name__=="__main__":
    main()