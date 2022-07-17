#!/usr/bin/python3

#----------------------------------------------------------------------------#
#
# Microbiome Insights | Snakemake pipeline for microbial characterization
#
#----------------------------------------------------------------------------#

# Module Imports
#----------------------------------------------------------------------------#

try: 
    import os, sys, re, itertools, time, argparse, copy,  configparser, subprocess, psutil, yaml
    import pandas as pd
    from pathlib import Path
except ModuleNotFoundError:
    print("")
    print("\033[93mLooks like a module was not found. Are you sure you have activated your conda environment?\033[0m")
    print("\033[93mSuggestion: `conda activate M5P`\033[0m")
    print("")
    sys.exit(1)

from datetime import datetime
now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

# This function looks for the Snakefile file in the path
def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), file)
    if not os.path.exists(sf):
        sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
        if not os.path.exists(sf):
            sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

# This function writes a template for the config file for the user is none is available
def get_template():
    m5P_config_template = {
    'bin_all': 'true',
    'database_dir': '~/M5P_databases',
    'experiment_type': 'metagenomics',
    'experimental_contrast': 'none',
    'experimental_design': '~treatment',
    'fastq_metagenomics': 'input_mg',
    'fastq_metatranscriptomics':' input_mt/',
    'genome': {
        'find_genome_relative': 'false',
        'genome_file': 'reference_file',
        'genome_url': 'reference_url',
        },
    'mag_annotation': 'dram',
    'max_memory': '1280000000000',
    'merged_prefix': 'MergedReads-001',
    'merged_reads': 'true',
    'metadata_path': 'metadata.txt',
    'threads': '2',
    'working_dir': 'tets./'
    }
    with open("template_M5P_config.yaml", 'w') as outfile:
        yaml.dump(m5P_config_template, outfile, default_flow_style=False)


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
    parser.add_argument("-i", "--fastq_metagenomics", type=str, help="Directory containing metagenomics fastq files (path)")
    parser.add_argument("-x", "--fastq_metatranscriptomics", type=str, help="Directory containing metatranscriptomics fastq files (path)")
#    parser.add_argument("-d", "--database_dir", type= str, help="Directory containing ATLAS databases (path)")
    parser.add_argument("-r", "--merged_reads", default= True, type=bool, help="Merge reads for co-assembly (True or False)")
    parser.add_argument("-m", "--metadata_path", type=str, help="Metadata file (path)")
    parser.add_argument("-t", "--threads", default= 2, type = int, help="The number of threads the pipeline is allowed to use (integer)")
    parser.add_argument("-M", "--max_memory", type = int, help = "The amount of memory provided to the assembler. Enter in byte format.")
    parser.add_argument("-e", "--experiment_type", type = str, default="metagenomics",  help = "Either metagenomics (default), metatranscriptomics, or both")
    parser.add_argument("-k", "--experimental_contrast", type = str, default="none",  help = "Experimental contrast as used in R formula")
    parser.add_argument("-g", "--experimental_design", type = str, default="~treatment",  help = "Experimental design as used in R formula")
    parser.add_argument("-p", "--merged_prefix", default="MergedReads-001", type=str, help="Prefix for merged reads files (string)")
    parser.add_argument("-b", "--bin_all", type=bool, default = True, help="Put all reads files in same bin group (bool)")
    parser.add_argument("-j", "--jobs", default = 2, type = str, help = "number of jobs")
    parser.add_argument("-c", "--configfile", default = "M5P_config.yaml", type=str, help = "optional yaml file containing all of the above configuration details.")
    parser.add_argument("-a", "--getconfig", action='store_true', help = "Prints a yaml template with the default configuration.")
    args = parser.parse_args()

    if args.getconfig == True:
        print("A template configuration (template_M5P_config.yaml), has been writen to the current directory")
        print("Please review and edit it before using it")
        get_template()
        sys.exit(1)

    if args.configfile:
        configfile = args.configfile

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
                    new_data["threads"] = val
                else:
                    raise Exception("Number of threads requested exceeds number of available cores.")
            elif arg == "max_memory":
                val = getattr(args, arg)
                if val <= psutil.virtual_memory()[1]:
                    max_memory = val
                    new_data["max_memory"] = val
                else:
                    raise Exception("Requested memory exceeds available memory, please change the value for the -M option")
            elif arg == "experiment_type":
                val = getattr(args, arg)
                if val in ["metagenomics", "metatranscriptomics", "both"]:
                    experiment_type = val
                    new_data["experiment_type"] = val
                else:
                    raise Exception("Experiment type not valid")
            elif arg == "working_dir":
                val = getattr(args, arg)
                new_data["working_dir"] = val
            elif arg == "experimental_contrast":
                val = getattr(args, arg)
                new_data["experimental_contrast"] = val
            elif arg == "experimental_design":
                val = getattr(args, arg)
                new_data["experimental_design"] = val
            elif arg == "fastq_metagenomics":
                val = getattr(args, arg)
                new_data["fastq_metagenomics"] = val
            elif arg == "fastq_metatranscriptomics":
                val = getattr(args, arg)
                new_data["fastq_metatranscriptomics"] = val
            elif arg == "merged_reads":
                val = getattr(args, arg)
                new_data["merged_reads"] = val
            elif arg == "merged_prefix":
                val = getattr(args, arg)
                new_data["merged_prefix"] = val
            elif arg == "metadata_path":
                val = getattr(args, arg)
                new_data["metadata_path"] = val
            elif arg == "bin_all":
                val = getattr(args, arg)
                new_data["bin_all"] = val
            elif arg == "jobs":
                val = getattr(args, arg)
                new_data["jobs"] = val

#            else: #Uncomment to allow 
#                new_data[arg] = val
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
    
    #Set paths (to pass on to snakemake as config)   
    snakepath = Path(get_snakefile())
    parent_dir = snakepath.parent.absolute()
 
    #Build snakemake command
    cmd = (
        "snakemake -s {snakefile} "
        "-j {jobs} --rerun-incomplete "
        "--use-conda "
        "--config parent_dir={parent_dir} "
        "--configfile {configfile}"
    ).format(
        snakefile=snakepath,
        jobs=args.jobs,
        parent_dir=parent_dir,
        configfile=configfile)
    print(f"Executing: {cmd}")

    #run snakemake command
    try:    
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        exit(1)

if __name__ == '__main__':
    main()
