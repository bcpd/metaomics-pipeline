#!/usr/bin/python3

# The shunning

import os
import gzip
import re
import sys
import argparse
import urllib.request
from urllib.error import HTTPError
import csv

def load_csv(filename):
    xopen = open
    if filename.endswith('.gz'):
        xopen = gzip.open
    with xopen(filename, "rt") as fp:
        r = csv.DictReader(fp)
        for row in r:
            yield row

def filterFile(inputFileName, outputFileName, unwanted_list):
    #input file reader
    infile = gzip.open(inputFileName, "rt")
    read = csv.reader(infile)
#    headers = read(read) # header

    #output file writer
    outfile = gzip.open(outputFileName, "wt")
    write = csv.writer(outfile)
#    write.writeheader() # write headers

    #for each row
    for row in read:
        #print(row)
        accession_id = row[6]
        accession_id = accession_id.split(' ')[0]
#        print(accession_id)
        if accession_id not in unwanted_list:
            # print(row)
            write.writerow(row)
    infile.close()
    outfile.close()

#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--genbank_cache', dest='genbank_dir',
        type=str, required=True,
        help='Genbank chache folder where genomes are initialy downloaded.')
    parser.add_argument('--grist_output_folder', dest='grist_dir', type=str, required=True,
        help='list of bracken files.')
    parser.add_argument('--output', dest='output',
        default='',required=False,
        help='Names of the genbank accessions that do not wor')
    args = parser.parse_args()
    genbank_dir = args.genbank_dir
    genbank_dir = genbank_dir.rstrip("\\/")
#    print(genbank_dir)
    grist_dir = args.grist_dir
    grist_dir = grist_dir.rstrip("\\/")
#   
    unavailable_genomes = []
    list_of_files = []
    for file in os.listdir(genbank_dir):
        if file.endswith('.info.csv'):
            list_of_files.append(file)
        else:
            continue
#    print(list_of_files)
    for input_file in list_of_files:
 #       print(input_file)
        rows = list(load_csv(genbank_dir + "/" + input_file))
        assert len(rows) == 1
        row = rows[0]
        ident = row['ident']
        url = row['genome_url']
        name = row['display_name']
        print(f"Reading {input_file} ")
        try:
            with urllib.request.urlopen(url) as response:
                content = response.read()
                print(f"Found a genome for {ident}")
        except urllib.request.HTTPError:
            print(f"Cannot download genome from URL:\n  {url}", file=sys.stderr)
            #raise Exception("Genbank genome not found")
            unavailable_genomes.append(ident)
    if len(unavailable_genomes) >0:
        print(f"\nThese genomes are not available: {unavailable_genomes}")

    list_of_prefetch = []
    for file in os.listdir(grist_dir + "/gather"):
        if file.endswith('.prefetch.csv.gz'):
            list_of_prefetch.append(file)
    #print(list_of_prefetch)
    for prefetch in list_of_prefetch:
        prefetch_in = grist_dir + "/gather/" + prefetch
        prefetch_out = grist_dir + "/temp_" + prefetch
        filterFile(prefetch_in, prefetch_out, unavailable_genomes)
        # replacement line
        replacement_command_1 = 'mv ' + prefetch_in + ' ' + grist_dir + '/original_' + prefetch
        print(replacement_command_1)
        # os.system(replacement_command_1)
        replacement_command_2 = 'mv ' + prefetch_out + ' ' + grist_dir + '/gather/' + prefetch
        print(replacement_command_2)
        # os.system(replacement_command_2)



# the main function 
main()

