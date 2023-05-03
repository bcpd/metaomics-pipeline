#!/usr/bin/python3

# The shunning

import os
import gzip
#import sys
import argparse
import urllib.request
from urllib.error import HTTPError
import csv


def does_genome_exist(accession):
    accsplit = accession.strip().split("_")
    assert len(accsplit) == 2, f"ERROR: '{accession}' should have precisely one underscore!"
    db, acc = accsplit
    if '.' in acc:
        number, version = acc.split(".")
    else:
        number, version = acc, '1'
    number = "/".join([number[p : p + 3] for p in range(0, len(number), 3)])
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/{number}"
    #print(f"opening directory: {url}", file=sys.stderr)
    with urllib.request.urlopen(url) as response:
        all_names = response.read()
        all_names = all_names.decode("utf-8")
        full_name = None
        for line in all_names.splitlines():
            if line.startswith(f'<a href='):
                name=line.split('"')[1][:-1]
                db_, acc_, *_ = name.split("_")
                if db_ == db and acc_.startswith(acc):
                    full_name = name
                    break
        if full_name is None:
            final_url = None
        else:
            final_url = "htt" + url[3:]
            final_url = (f"{url}/{full_name}/{full_name}_genomic.fna.gz")
        try:
            with urllib.request.urlopen(final_url) as response:
                content = response.read()
                return(True)
        except urllib.request.HTTPError:
            #print(f"Cannot download genome from URL:\n  {final_url}", file=sys.stderr)
            return(False)


#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--grist_output_folder', dest='grist_dir', type=str, required=True, help='Grist output folder')
    args = parser.parse_args()
    grist_dir = args.grist_dir
    grist_dir = grist_dir.rstrip("\\/")
#   
    unavailable_genomes = []

    list_of_prefetch = []
    for file in os.listdir(grist_dir + "/gather"):
        if file.endswith('.prefetch.csv.gz'):
            list_of_prefetch.append(file)

    for prefetch in list_of_prefetch:
        print(f"Working on file {prefetch}")
        prefetch_in = grist_dir + "/gather/" + prefetch
        prefetch_out = grist_dir + "/temp_" + prefetch
        infile = gzip.open(prefetch_in, "rt")
        read = csv.reader(infile)
        outfile = gzip.open(prefetch_out, "wt")
        write = csv.writer(outfile)
        for row in read:
            ident = row[6]
            accession_id = ident.split(' ')[0]
            #print(accession_id)
            if accession_id == "match_name" or does_genome_exist(accession_id):
#                print(f"genome {accession_id} exist")
                write.writerow(row)
            else:
                print(f"A genome for {accession_id} is not available")
                if accession_id not in unavailable_genomes:
                    unavailable_genomes.append(accession_id)
                else:
                    continue
        infile.close()
        outfile.close()
        # replacement line
        replacement_command_1 = 'mv ' + prefetch_in + ' ' + grist_dir + '/original_' + prefetch
        #print(replacement_command_1)
        os.system(replacement_command_1)
        replacement_command_2 = 'mv ' + prefetch_out + ' ' + grist_dir + '/gather/' + prefetch
        #print(replacement_command_2)
        os.system(replacement_command_2)

    if len(unavailable_genomes) >0:
        print(f"\nThese genomes are not available: {unavailable_genomes}")


# the main function 
main()