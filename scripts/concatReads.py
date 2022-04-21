#!/usr/bin/env python

#---------------------------------------------------------------------------
# concatReads.py
# Usage: run concatReads.py to concatonate reads 
# Supply a directory and the script will find and 
# Merge gzipped read files (fastq.gz or fq.gz)
# Optionally, will downsample using bbmap (must be installed)
# Author: David Levy-Booth
#---------------------------------------------------------------------------

import os, sys, argparse, subprocess


#Arguments and setting
#---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True,
                    help="Input directory containing read files")
parser.add_argument("-o", "--output", type=str, default="allReads",
                    help="Output file basename (optional)")
parser.add_argument("-d", "--downrate", type=float, default=0,
                    help="Downsampling rate (optional)")
parser.add_argument("-r", "--downreads", type=int, default=0,
                    help="Downsampling read depth (optional)")
args = parser.parse_args()

#Settings
rFormat= {"fwd": "R1", "rev": "R2"}
qFormat= (".fastq.gz", ".fq.gz")


#Functions
#---------------------------------------------------------------------------


def collectReads(inputdir, rFormat, qFormat):
    #Extract files from directory
    #Split on rFormat strings for 
    #fwd and rev reads. Return as lists
    if not os.path.isdir(inputdir):
        print("ERROR: --input not a directory")
    else:
        reads = [f for f in os.listdir(inputdir) if f.endswith(qFormat)]
        reads.sort()
        #Filter forward and reverse reads
        f = [str(os.path.join(inputdir, f)) for f in reads if rFormat["fwd"] in f]
        r = [str(os.path.join(inputdir, r)) for r in reads if rFormat["rev"] in r]
        return f, r

def buildOut(rd, inputdir, output, rFormat, fqSuffix):
    outRead = str(os.path.join(inputdir, f"{output}_{rFormat[rd]}{fqSuffix}"))
    return outRead


def concatReads(fwd, rev, rFormat, fqSuffix):
    #build output files
    outFwd = buildOut("fwd", args.input, args.output, rFormat, fqSuffix)
    outRev = buildOut("rev", args.input, args.output, rFormat, fqSuffix)
    #Call concatonate (cat)
    if os.path.isfile(outFwd):
        print(f"ERROR: output file {outFwd} already exists")
    else:
        fwdCall = ["cat"] + fwd + [">"] + [f"{outFwd}"]
        revCall = ["cat"] + rev + [">"] + [f"{outRev}"]

        os.system(" ".join(fwdCall))
        os.system(" ".join(revCall))

    return outFwd, outRev


def downsampleReads(fwd, rev, rFormat, fqSuffix):
    seed=2022
    #build 'output' files
    outFwd = buildOut("fwd", args.input, args.output, rFormat, fqSuffix)
    outRev = buildOut("rev", args.input, args.output, rFormat, fqSuffix)
    outFwdSub = buildOut("fwd", args.input, args.output + f"Sub{args.downrate}", rFormat, fqSuffix)
    outRevSub = buildOut("rev", args.input, args.output + f"Sub{args.downrate}", rFormat, fqSuffix)

    #Create process to cat input files
    #input process to input of reformat.sh
    if os.path.isfile(outFwdSub):
        print(f"ERROR: output file {outFwdSub} already exists")
    #Forward Read 
    fps = subprocess.run(["cat"] + fwd, check=True, capture_output=True)
    subprocess.run(["reformat.sh", 
    f"in=stdin.fq.gz",
    f"out1={outFwdSub}",
    f"samplerate={args.downrate}",
    "int=false",
    "overwrite=true"],
    input=fps.stdout)

    #Reverse Read 
    fps = subprocess.run(["cat"] + rev, check=True, capture_output=True)
    subprocess.run(["reformat.sh", 
    f"in=stdin.fq.gz",
    f"out1={outRevSub}",
    f"samplerate={args.downrate}",
    "int=false",
    "overwrite=true",
    f"sampleseed={seed}"],
    input=fps.stdout)

    # #Forward Read 
    # subprocess.run(["reformat.sh", 
    # f"in1={outFwd}",
    # f"in2={outRev}",
    # f"out1={outFwdSub}",
    # f"out2={outFwdSub}",
    # f"samplerate={args.downrate}",
    # "int=false",
    # "overwrite=true"])


    #repair.sh

#Main implementation
#---------------------------------------------------------------------------
if __name__ == '__main__':
    #Create seperate lists for "R1" and "R2" files
    fwd, rev = collectReads(args.input, rFormat, qFormat)

    #downsample
    if args.downrate:
        downsampleReads(fwd, rev, rFormat, ".fastq.gz")
    else:
        #call the bash 'cat' function
        outFwd, outRev = concatReads(fwd, rev, rFormat, ".fastq.gz")
