#!/usr/bin/env python

# -*- coding: utf-8 -*-

import os
import sys
import time
import argparse
import yaml
from datetime import datetime


usage = '''
Usage:
     pyton create_grist_config_file.py
         -s <samples folder>
         -d <grist database folder>
         -o <grist_configuration_file.yaml>
         -w <working directory>
         -m <maximum_memory
'''


def main():
    # Import sample names using pandas
    parser = argparse.ArgumentParser(description='Grist configuration setup')
    parser.add_argument("-s", dest="samples_folder", type=str,
                        help="Folder with metatranscriptomics fastq(path)")
    parser.add_argument("-o", dest="output_fp", type=str,
                        default="grist_config.yaml",
                        help="Directory containing fastq files (path)")
    parser.add_argument("-w", dest="working_dir", type=str,
                        help="Workinf directory")
    parser.add_argument("-m", dest="max_memory", type=int, default=10e9,
                        help="Maximum RAM memory")

    parser.add_argument("-d", dest="databases_folder",
                        default='~/databases/grist/', type=str,
                        help="Directory containing grist databases (path)")
    (opts, args) = parser.parse_args()

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print(now)
    print('Initializing...')
    print('Using files in: ' + opts.samples_folder)
    print('Databases in: ' + opts.samples_folder)
    print('Writing to: ' + opts.output_fp)

    # Find samples and modify them if needed
    sample_list = os.listdir(opts.samples_folder)
    if len(sample_list) == 0:
        print("There were no files in the samples folder ")
        print(usage)
        sys.exit(0)

    samples = {}

    # Create temporary folder
    working_dir = opts.working_dir
    working_dir = working_dir.rstrip('/', 1)
    grist_dir = working_dir + '/grist'
    raw_files_dir = working_dir + '/grist/raw'

    os.system('mkdir -p ' + working_dir + '/grist')
    os.system('mkdir -p ' + working_dir + '/grist/raw')

    for i in sample_list:
        if i.endswith('.fastq.gz'):
            if i.endswith('_1.fastq.gz'):
                original_basename = i.rsplit('_', 1)[0]
                file_basename = original_basename.replace(".", "_")
                file_basename = file_basename.replace("-", "_")
                samples.append(opts.sample_folder + '/' + file_basename)
                os.system('cp ' +
                          opts.sample_folder + '/' +
                          original_basename + '_1.fastq.gz ' +
                          raw_files_dir + '/' +
                          file_basename + '_1.fastq.gz')
                os.system('cp ' +
                          opts.sample_folder + '/' +
                          original_basename + '_2.fastq.gz ' +
                          raw_files_dir + '/' +
                          file_basename + '_2.fastq.gz')
            elif i.endswith('_R1_001.fastq.gz'):
                original_basename = i.rsplit('_', 2)[0]
                file_basename = original_basename.replace(".", "_")
                file_basename = file_basename.replace("-", "_")
                samples.append(opts.sample_folder + '/' + file_basename)
                os.system('cp ' +
                          opts.sample_folder + '/' +
                          original_basename + '_R1_001.fastq.gz ' +
                          raw_files_dir + '/' +
                          file_basename + '_1.fastq.gz')
                os.system('cp ' +
                          opts.sample_folder + '/' +
                          original_basename + '_R2_001.fastq.gz ' +
                          raw_files_dir + '/' +
                          file_basename + '_2.fastq.gz')
            else:
                print("There were no fastq files in the samples folder ")
                print("I expected fastq.gz files")
                print(usage)
                sys.exit(0)
        else:
            continue
    print(samples)

    # Find databases in folder
    # If no databases are found report an error

    databases_folder = opts.databases_folder
    databases_folder = databases_folder.rstrip('/', 1)
    dbs_list = os.listdir(databases_folder)
    # Check that the folder is not empty
    if len(dbs_list) == 0:
        print("There were no databases in the folder ")
        print(usage)
        sys.exit(0)
    # Create empty list to store names of databases
    my_grist_databases = []
    for db in dbs_list:
        if db.endswith('.zip'):
            my_grist_databases.append(opts.databases_folder + '/' + db)
            continue
        else:
            print("I was expecting zip files for databases files")
            print(usage)
            sys.exit(0)

    # Write results into yaml file
    # Expected grist config_file
    """
    # example config file #
    samples:
    - Mock2_T0_1_S49
    - Mock2_T0_2_S50
    - Mock2_T0_3_S51

    outdir: outputs.tutorial.mtx3/

    metagenome_trim_memory: 10e9

    sourmash_databases:
    - db/gtdb-rs202.genomic-reps.k31.sbt.zip

    prevent_sra_download: true
    """

    grist_data = {
            'samples': samples,
            'outdir': opts.output_fp,
            'sourmash_databases': my_grist_databases,
            'prevent_sra_download': 'true', 
            'outdir': grist_dir
            }

    with open(opts.output_fp, 'w') as outfile:
        yaml.dump(grist_data, outfile, default_flow_style=False)

if __name__ == ' __main__':
    main()
