#!/usr/bin/python3

#----------------------------------------------------------------------------#
#
# Microbiome Insights | Snakemake pipeline for microbial characterization
#
#----------------------------------------------------------------------------#

def COLLECT_ALL_INPUT():
    INPUTS = []
    if config["experiment_type"] == "metagenomics":
        #Config outputs for metagenomics
        INPUTS.append(os.path.join(working_dir, "samples.tsv")) #rule atlas_init
        if merged_reads:
            INPUTS.append(os.path.join(working_dir, "logs/concatReads.log"))
        INPUTS.append(os.path.join(working_dir, "logs/formatSamples.log"))
        INPUTS.append(os.path.join(working_dir, "finished_genecatalog"))
        INPUTS.append(os.path.join(working_dir, "finished_genomes"))
        INPUTS.append(os.path.join(working_dir, "logs/DRAM_copy_results.log")) # DRAM
        INPUTS.append(os.path.join(working_dir, "logs/run_bakta.log")) # bakta
        INPUTS.append(os.path.join(working_dir, "logs/run_referenceseeker.log")) # referenceseeker
        INPUTS.append(os.path.join(working_dir, "logs/Atlas_metagenomics_cleanup.log")) # referenceseeker
    if config["experiment_type"] == "metatranscriptomics":
        #Config outputs for metatranscriptomics
        INPUTS.append(os.path.join(working_dir, "grist/reports")) #grist
        INPUTS.append(os.path.join(working_dir, "salmon/DESeq2.tsv") # praxis

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

def COLLECT_SALMON_REFERENCE():
    INPUTS=[]
    if "metagenomics" in config["experiment_type"]:
        #Config outputs for metagenomics
        INPUTS.append(os.path.join(working_dir, "Genecatalog/gene_catalog.fna")) #rule atlas_genecatalog
    elif config["experiment_type"] == "metatranscriptomics":
        #Config outputs for metatranscriptomics
        INPUTS.append(os.path.join(working_dir, "grist/reports/*.fna")) #grist
