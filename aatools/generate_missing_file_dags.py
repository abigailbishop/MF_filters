from collections.abc import Iterable

def generate_dag(runlist: Iterable, save_dir: str, station: int) -> None: 
    """
    Given the runs in the `runlist` and the `station` to be analyzed, save 
      dagman commands to a file located in `save_dir`. 

    Parameters
    ----------
    runlist : Iterable
        List/tuple/etc of runs to create a dagman for.
    save_dir : str
        Path to the directory where you'd like to save the new dag file to.
    station : int
        The integer value corresponding to the station you are analyzing. 

    Returns
    -------
    None
    """

    # Open the file the new dagman commands will be saved to
    print("Saving new dag to:",save_dir+f"A{station}_missing.dag")
    file = open(save_dir+f"A{station}_missing.dag", "w")

    # Add a dagman command to the `file` for each run in the `runlist`
    for run in runlist: 
        file.write(f"JOB job_ARA_S2_R{run} ARA_job.sub\n")
        file.write(f'VARS job_ARA_S2_R{run} station="{station}" run="{run}"\n\n')

    # Close the `file`
    file.close()

    return


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

    # Collect all runs from `args.runslit` that are missing from 
    #   `args.output_dir` and are between `args.first_run` and `args.last_run`.
    from find_missing_outputs import find_missing_outputs
    missing_runs = find_missing_outputs(
        args.runlist, args.output_dir, args.first_run, args.last_run)
    
    # Create the file with dag submission commands for each missing run. 
    generate_dag(missing_runs, args.save_dir, args.station)

    return

if __name__=="__main__":
    main()