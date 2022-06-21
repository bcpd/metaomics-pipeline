#!/usr/bin/python3

#----------------------------------------------------------------------------#
#
# Microbiome Insights | Snakemake pipeline for microbial characterization
#
#----------------------------------------------------------------------------#

def COLLECT_ALL_INPUT():
    INPUTS = []
    INPUTS.append(os.path.join(working_dir, "samples.tsv")) #rule atlas_init
    if merged_reads:
        INPUTS.append(os.path.join(working_dir, "logs/concatReads.log"))
    INPUTS.append(os.path.join(working_dir, "logs/formatSamples.log"))
    INPUTS.append(os.path.join(working_dir, "finished_genecatalog"))
    INPUTS.append(os.path.join(working_dir, "finished_genomes"))

    return INPUTS

def COLLECT_INIT_INPUT():
    INPUTS = []
    INPUTS.append(fastq_metagenomics)
    if merged_reads:
        INPUTS.append(os.path.join(working_dir, "logs/concatReads.log"))
    return INPUTS

def COLLECT_FORMAT_ARGS():
    if bin_all:
        binarg = "--all"
    else:
        binarg = ""
    if metadata_path:
        return f"-s {os.path.join(working_dir, 'samples.tsv')} -m {metadata_path} {binarg}"
    else:
        return f"-s {os.path.join(working_dir, 'samples.tsv')} {binarg}"
