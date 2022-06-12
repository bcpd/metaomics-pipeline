#!/usr/bin/python3

#----------------------------------------------------------------------------#
#
# Microbiome Insights | Snakemake pipeline for microbial characterization
#
#----------------------------------------------------------------------------#


# Module Imports
#----------------------------------------------------------------------------#
import M5P

import os, sys, re, itertools, time, argparse, psutil, copy, yaml, configparser, subprocess
from snakemake.utils import validate

import pandas as pd
from pathlib import Path
from datetime import datetime
now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

# This function looks for the Snakefile file in the path
def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf



def main():
    print("\033[92m                                             ,,,.       ,,,                 \033[0m")
    print("\033[92m                                         ,      %%%%     ,,,   ,            \033[0m")
    print("\033[92m                                      ,,      %%%%%       */**    ,,        \033[0m")
    print("\033[92m                                   .,         %%%%                   ,      \033[0m")
    print("\033[92m        __  __ _            _    _                 (%&&&&&&%%          ,    \033[0m")
    print("\033[92m       |  \\/  (_)__ _ _ ___| |__(_)___ _ __  ___   %%&&&&&&%  %%%%%     ,,  \033[0m")
    print("\033[92m       | |\\/| | / _| '_/ _ \\ '_ \\ / _ \\ '  \\/ -_)              %%%%%     ,,\033[0m")
    print("\033[92m       |_|  |_|_\\__|_| \\___/_.__/_\\___/_|_|_\\___|                          ,\033[0m")
    print("\033[92m ___                 _           _       _                                  ,\033[0m")
    print("\033[92m|_ _|  _ __    ___  (_)   __ _  | |__   | |_   ___          %%%%            ,\033[0m")
    print("\033[92m | |  | '_ \\  / __| | |  / _` | | '_ \\  | __| / __|          %%%%           ,\033[0m")
    print("\033[92m | |  | | | | \\__ \\ | | | (_| | | | | | | |_  \\__ \\             /%         ,\033[0m")
    print("\033[92m|___| |_| |_| |___/ |_|  \\__, | |_| |_|  \\__| |___/                        ,\033[0m")
    print("\033[92m                          |___/        **,         %%%%%%        %%%%%   ,, \033[0m")
    print("\033[92m                                ,,   ,*/*             /%        &&&&&&  ,,  \033[0m")
    print("\033[92m   MBI Microbial                  ,   ,    %%%&&&              %&&&&&  ,    \033[0m")
    print("\033[92m   Mixtures Metagenomics and       .,      %%&&@&&&%    ***,  %%%%%  ,      \033[0m")
    print("\033[92m   Metatranscriptomics pipeline v1.0           &&&%%    **/*,     ,,        \033[0m")
    print("\033[92m                                         ,,          ,,,       ,,           \033[0m")
    print("\033[92m   (Perfect, no notes)                       .   ,      ,,,,                \033[0m")
    time.sleep(1)


    parser = argparse.ArgumentParser(description='pipeline configuration')
    parser.add_argument("-w", "--working_dir", default="./", type=str, help ="working directory (path)")
    parser.add_argument("-i", "--fastq_dir", type=str, help="Directory containing fastq files (path)")
    parser.add_argument("-d", "--database_dir", type= str, help="Directory containing ATLAS databases (path)")
    parser.add_argument("-r", "--merged_reads", default= True, type=bool, help="Merge reads for co-assembly (True or False)")
    parser.add_argument("-m", "--metadata_path", type=str, help="Metadata file (path)")
    parser.add_argument("-t", "--threads", default= 2, type = int, help="The number of threads the pipeline is allowed to use (integer)")
    parser.add_argument("-M", "--max_memory", type = int, help = "The amount of memory provided to the assembler. Enter in byte format.")
    parser.add_argument("-p", "--merged_prefix", default="MergedReads-001", type=str, help="Prefix for merged reads files (string)")
    parser.add_argument("-b", "--bin_all", type=bool, default = True, help="Put all reads files in same bin group (bool)")
    parser.add_argument("-j", "--jobs", default = 2, type = str, help = "number of jobs")
    parser.add_argument("-c", "--configfile", default = "M5P_config.yaml", type=str, help = "optional yaml file containing all of the above configuration details.")
    args = parser.parse_args()

    if args.configfile:
        configfile = args.configfile
        try:
            validate(config, configfile)
        except:
            raise ImportError(f"ERROR: confirm config file {configfile} path and formatting")


    # Open and create an image of the Snakemake config
    stream = open(configfile, "r")
    original_data = yaml.load(stream, yaml.FullLoader)
    new_data = copy.deepcopy(original_data)

    # Check for corect values for the arguments
    for arg in vars(args):
        if getattr(args, arg) is not None:
            if arg == "threads":
                val = getattr(args, arg)
                if val <= psutil.cpu_count():
                    threads = val
                else:
                    raise Exception("Number of threads requested exceeds number of available cores.")
            elif arg == "max_memory":
                val = getattr(args, arg)
                if val <= psutil.virtual_memory()[1]:
                    max_memory = val
                else:
                    raise Exception("Requested memory exceeds available memory")
#           elif getattr(args, arg) is not None:
#                if arg == "a flag":
#                    val = getattr(args, arg)
#                if test_for_that flag:
#                    flag = val
#                else:
#                    raise Exception("")


    # Rewrite config file if changes were made to image
    if original_data != new_data:
        with open(configfile, 'w') as yaml_file:
            yaml_file.write(yaml.dump(new_data, default_flow_style=False))

    #Build snakemake command
    cmd = (
        "snakemake -s {snakefile} "
        "-j {jobs} --rerun-incomplete "
        "--use-conda "
        "--configfile {configfile}"
    ).format(
        snakefile=get_snakefile(),
        jobs=jobs,
        configfile=configfile),
    logging.info(f"Executing: {cmd}")

    #run snakemake command
    try:    
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        exit(1)

if __name__ == ' __main__':
    main()