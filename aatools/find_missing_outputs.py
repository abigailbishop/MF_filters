from collections.abc import Iterable

def get_runs_from_runlist(
    list_name: str, first_run: int = 0, last_run: int = 100000
) -> tuple :
    """
    Open the runlist at the provided location `list_name` and extract all runs.

    Parameters
    ----------
    list_name : str
        Path to the run list.
    first_run : int
        The first run you'd like to be pulled from the run list. 
        If `first_run` is greater than the maximum run in the run list, 
          no runs should be returned.
    last_run : int
        The last run in the run list you'd like to be extracted from the list.
        If `last_run` is less than the minimum run in the run list, 
          no runs should be returned.

    Returns
    -------
    runs : tuple
        A tuple containing every run from the run list (between `first_run` and
          `last_run`) as an integer
    """

    # Open the run list 
    list_file = open(list_name, "r")

    # Collect all runs in the run list
    runs = []
    for lines in list_file:
        line = lines.split()
        run_num = int(line[0])
        if (run_num > first_run) and (run_num < last_run): 
            runs.append(run_num)
        del line, run_num

    # Close the run list
    list_file.close()

    # Return all runs as a tuple of integers
    return tuple(runs)

def get_runs_from_output_dir(output_dir: str) -> tuple:
    """
    Find and report all the runs located in the `output_dir`

    Parameters
    ----------
    output_dir : str
        The path to the output directory you'd like to find all runs in.

    Returns
    -------
    runs : tuple
        A tuple containing every run located in the output directory as an int.
    """

    # Collect all files in the output directory
    from os import listdir
    output_files = listdir(output_dir)
    
    # Extract the run associated with each file in the output directory
    runs = []
    for file in output_files: 
        # sub_info_full_A2_R18620.h5
        run_num = int(file.split(".")[-2].split("R")[-1])
        runs.append(run_num)

    # Return all the identified runs as a tuple of integers
    return tuple(runs)

def find_missing_files(input_runs: Iterable, output_runs: Iterable) -> tuple:
    """
    Identify all runs in `output_runs` that are not in `input_runs`.

    Parameters
    ----------
    input_runs : Iterable
        List/tuple/etc of all expected runs from the run_list.
    output_runs : Iterable
        List/tuple/etc of all runs that have an output file. 

    Outputs
    -------
    missing_runs : tuple
        Tuple of all runs in `input_runs` that are missing from `output_runs.
    """

    # For each run in `input_runs`, check that it is in `output_runs`. If not,
    #   log that run in `missing_runs.`
    missing_runs = []
    output_runs = set(output_runs)
    for input_run in input_runs: 
        if input_run not in output_runs: 
            missing_runs.append(input_run)

    # Return missing_runs as an immutable tuple
    return tuple(missing_runs)

def find_missing_outputs(
    runlist: str, output_dir: str, first_run: int, last_run: int
) -> tuple:
    """
    Find all runs from the `runlist` (between `first_run` and `last_run`) that 
      are missing from the `output_dir`.

    Parameters
    ----------
    runlist : str
        Path to the list of all runs.
    output_dir : str
        Path to the directory containing output files for each run.
    first_run : int
        The first run in the `runlist` intended for analysis. 
    last_run : int
        The last run in the `runlist` intended for analysis. 

    Returns
    -------
    missing_files : tuple
        All runs for which an associated file was not found in `output_dir` 
          even though the run exists in the `runlist` and is between 
          `first_run` and `last_run`. 
    """

    # Get all input runs from the `runlist` between `first_run` and `last_run`
    input_runs = get_runs_from_runlist(
        runlist, first_run, last_run)
    
    # Get all output runs from the `output_dir`
    output_runs = get_runs_from_output_dir(output_dir)

    # Identify all runs from the `runlist` that are missing from the `output_dir`
    missing_files = sorted(find_missing_files(input_runs, output_runs))
    
    # Share statistics with the user
    print("Number of input runs:", len(input_runs))
    print("Number of output runs:", len(output_runs))
    print("Number of missing files:", len(missing_files))
    print()
    print(missing_files)

    # Return all runs from the `runlist` without files in the `output_dir`
    return tuple(missing_files)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runlist', type=str,
        help="Absolute path to runlist")
    parser.add_argument('output_dir', type=str,
        help="Absolute path to output files.")
    parser.add_argument('--first_run', type=int, default = 0,
        help="First run being analyzed")
    parser.add_argument('--last_run', type=int, default = 100000,
        help="Last run being analyzed")
    args = parser.parse_args()

    find_missing_outputs(
        args.runlist, args.output_dir, 
        first_run=args.first_run, last_run=args.last_run,)