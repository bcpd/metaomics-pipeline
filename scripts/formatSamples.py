#!/usr/bin/env python

#---------------------------------------------------------------------------
# formatSamples.py
# Usage: run formatSamples.py to modify the
# samples.tsv file to merge bin groups according to 
# an optional metadata.tsv file
# 
# Author: David Levy-Booth
#---------------------------------------------------------------------------

import os, sys, argparse, subprocess
import yaml
import pandas as pd


#Arguments and setting
#---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--samples", type=str, required=True,
                    help="Input samples.tsv file path")
# parser.add_argument("-c", "--config", type=str, required=True,
#                     help="Input config.yaml file path")
parser.add_argument("-m", "--metadata", type=str,
                    help="Metadata file (optional)")
parser.add_argument("-a", "--all", action='store_true',
                    help="Put all samples in same bin group")
args = parser.parse_args()

#verifypaired= "false"


#Main implementation
#---------------------------------------------------------------------------
if __name__ == '__main__':
    #Load samples
    samples = pd.read_csv(args.samples, 
                          sep="\t")

    #If metadata file, use the SampleName column to 
    #match the rows against the samples.tsv file
    #then change BinGroups to match the metadata file
    if args.metadata:
        metadata = pd.read_csv(args.metadata, 
                          sep="\t")
        
        metadata = metadata.set_index('SampleName')
        metadata = metadata.reindex(index=samples['Full_Name'])
        metadata = metadata.reset_index()

        samples["BinGroup"] = metadata["BinGroup"]

    #As an alterative strategy, set only one bin group 
    if args.all:
        samples["BinGroup"] == "BinGroup1"

    #Overwrite file
    samples.to_csv(args.samples, sep="\t", index = False)

    # #Read yaml, add verification argument
    # with open(args.config, 'r') as file:
    #     config = yaml.safe_load(file)
    #     config["verifypaired"] = verifypaired
    # #Overwrite yaml
    # with open(args.config, 'w') as file:
    #     yaml.dump(config, 
    #             file, 
    #             default_style=None, 
    #             default_flow_style=False, 
    #             sort_keys=False)