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

def url_for_accession(accession):
    accsplit = accession.strip().split("_")
    assert len(accsplit) == 2, f"ERROR: '{accession}' should have precisely one underscore!"

    db, acc = accsplit
    if '.' in acc:
        number, version = acc.split(".")
    else:
        number, version = acc, '1'
    number = "/".join([number[p : p + 3] for p in range(0, len(number), 3)])
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/{number}"
    print(f"opening directory: {url}", file=sys.stderr)

    with urllib.request.urlopen(url) as response:
        all_names = response.read()

    print("done!", file=sys.stderr)

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
        return None
    else:
        url = "htt" + url[3:]
        return (
            f"{url}/{full_name}/{full_name}_genomic.fna.gz")


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
    parser.add_argument('--grist_output_folder', dest='grist_dir', type=str, required=True, help='list of bracken files.')
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
        rows = list(load_csv(grist_dir + "/gather/" + prefetch))
        for row in rows:
            ident = row['match_name']
            accession = ident.split(' ')[0]
            accsplit = accession.strip().split("_")
            assert len(accsplit) == 2, f"ERROR: '{accession}' should have precisely one underscore!"
            db, acc = accsplit
            if '.' in acc:
                number, version = acc.split(".")
            else:
                number, version = acc, '1'
            number = "/".join([number[p : p + 3] for p in range(0, len(number), 3)])
            url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/{number}"
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

            #print(f"opening directory: {url}", file=sys.stderr)
            try:
                with urllib.request.urlopen(final_url) as response:
                    content = response.read()
                    print(f"Found a genome for {ident}")
            except urllib.request.HTTPError:
                print(f"Cannot download genome from URL:\n  {final_url}", file=sys.stderr)
                #raise Exception("Genbank genome not found")
                if ident not in unavailable_genomes:
                    unavailable_genomes.append(ident)

    if len(unavailable_genomes) >0:
        print(f"\nThese genomes are not available: {unavailable_genomes}")

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

    if len(unavailable_genomes) >0:
        print(f"\nThese genomes are not available: {unavailable_genomes}")




# the main function 
main()

