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
import pandas as pd
from pathlib import Path
from datetime import datetime
now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

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
    parser.add_argument("-w", "--working_dir", default = "./", type = str, help = "working directory (path)")
    parser.add_argument("-i", "--fastq_dir", type = str, help = "Directory containing fastq files (path)")
    parser.add_argument("-d", "--database_dir", type = str, help = "Directory containing ATLAS databases (path)")
    parser.add_argument("-r", "--merged_reads", default = True, type = bool, help = "Merge reads for co-assembly (True or False)")
    parser.add_argument("-m", "--metadata_path", type = str, help = "Metadata file (path)")
    parser.add_argument("-t", "--threads", default = 2, type = int, help = "The number of threads the pipeline is allowed to use (integer)")
    parser.add_argument("-p", "--merged_prefix", default = "MergedReads-001", type = str, help = "Prefix for merged reads files (string)")
    parser.add_argument("-b", "--bin_all", type = bool, default = True, help = "Put all reads files in same bin group (bool)")
    parser.add_argument("-c", "--configfile", type = str, help = "optional yaml file containing all of the above configuration details.")
    args = parser.parse_args()

    if args.configfile:
        configfile = args.configfile

    threads = 16
    snakefile = "~/DAVIDDEV/metaomics-pipeline/M5P/Snakefile"
    #Actually running snakemake
    cmd = f"snakemake -s {snakefile} -j{threads} --rerun-incomplete"
    # snakemake ~/DAVIDDEV/metaomics-pipeline/M5P/Snakefile -j16 --rerun-incomplete
    try:    
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        exit(1)


