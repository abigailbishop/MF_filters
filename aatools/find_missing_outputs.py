def get_runs_from_runlist(list_name, first_run = 0, last_run=100000):

    list_file =  open(list_name, "r")
    runs = []
    for lines in list_file:
        line = lines.split()
        run_num = int(line[0])
        if (run_num > first_run) and (run_num < last_run): 
            runs.append(run_num)
        del line, run_num
    list_file.close()

    return tuple(runs)

def get_runs_from_output_dir(output_dir):
    from os import listdir
    
    runs = []
    output_files = listdir(output_dir)
    for file in output_files: 
        # sub_info_full_A2_R18620.h5
        run_num = int(file.split(".")[-2].split("R")[-1])
        runs.append(run_num)

    return tuple(runs)

def find_missing_files(input_runs, output_runs):

    missing_runs = []
    output_runs = set(output_runs)
    for input_run in input_runs: 
        if input_run not in output_runs: 
            missing_runs.append(input_run)

    return tuple(missing_runs)

def find_missing_outputs(runlist, output_dir, first_run, last_run):

    input_runs = get_runs_from_runlist(
        runlist, first_run, last_run)
    output_runs = get_runs_from_output_dir(output_dir)
    missing_files = sorted(find_missing_files(input_runs, output_runs))
    print(len(input_runs), len(output_runs), len(missing_files))
    print(missing_files)

    return missing_files

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