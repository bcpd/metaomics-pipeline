#!/usr/bin/python3

#----------------------------------------------------------------------------#
#
# Microbiome Insights | Snakemake pipeline for microbial characterization
#
#----------------------------------------------------------------------------#


# Module Imports
#----------------------------------------------------------------------------#

import os, sys, re, itertools, argparse, psutil, copy, yaml, configparser
import pandas as pd
from pathlib import Path
from datetime import datetime
now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")


# Configuration 
#----------------------------------------------------------------------------#

# EXPERIMENT_TYPE = config["EXPERIMENT_TYPE"]

# MAG_ANNOTATION = config["MAG_ANNOTATION"]

# ASSEMBLY_STRATEGY = config["ASSEMBLY_STRATEGY"]

# GENOME_REFERENCES = config["GENOME_REFERENCES"]

# ANNOTATION_DATABASES = config["ANNOTATION_DATABASES"]

# # Experimental design


# samples_table = pd.read_csv(config["samples"], sep="\t")
# sampleID_list = samples_table["SampleID"]


# EXP_METADATA = config["METADATA"]
# EXP_DESIGN = config["EXP_DESIGN"]
# EXP_CONTRAST = config["EXP_CONTRAST"]

# # Computational resources
# N_CORES = config["N_CORES"]
# MAX_MEMORY = config["MAX_MEMORY"]

# #

# METAGENOMICS_INPUT_FOLDER = config["MAX_MEMORY"]
# METATRANSCRIPTOMICS_INPUT_FOLDER = config["MAX_MEMORY"]

# OUTPUT_FOLDER = config["OUTPUT_FOLDER"]



# Begin Snakemake 

# # ---- Targets rules
# include: "workflow/rules/targets/targets.smk"

# # ---- Quality control rules
# include: "workflow/rules/qc/atlas_init.smk"
# include: "workflow/rules/qc/atlas_qc.smk"

# # ---- Assembly rules
# include: "workflow/rules/assembly/atlas_assembly.smk"
# include: "workflow/rules/assembly/atlas_binning.smk"
# include: "workflow/rules/assembly/atlas_coassembly.smk"

# # ---- Contig annotation rules
# include: "workflow/rules/annotation/atlas_annotation.smk"
# include: "workflow/rules/annotation/dram_annotation.smk"
# include: "workflow/rules/annotation/bakta_annotations.smk"
# include: "workflow/rules/annotation/atlas_taxonomy.smk"
# #include: "workflow/rules/annotation/integrate_annotations.smk"


# # ---- Metatranscriptomics rules
  
# include: "workflow/rules/metatranscriptomics/praxis_quality_control.smk"

#  # Use this option is a link to an NCBI genome is provided
# if GENOME_REFERENCES == "referemnce_url":
#     include: "workflow/rules/metatranscriptomics/praxis_download.smk"
#     annotate: false
# elif GENOME_REFERENCES == "referemnce_file":     
#       annotate: false
# # Use this option if you want to make a denove metranscriptome assembly
# elif GENOME_REFERENCES == "assemble_metatranscriptome":  
#     include: "workflow/rules/metatranscriptomics/praxis_assemble.smk"
# # Use this option if you want to find close relatives with grist
# else GENOME_REFERENCES == "find_relatives":
#     include: "workflow/rules/metatranscriptomics/grist.smk"
#      annotate: false

# include: "workflow/rules/metatranscriptomics/praxis_index_align.smk"
# include: "workflow/rules/metatranscriptomics/praxis_quantify.smk"
# include: "workflow/rules/metatranscriptomics/praxis_transcript_annotate.smk"

#M5P CONFIG
fastq_dir     = "/home/mixtures/test_data_d1/mtg"
database_dir  = "/home/mixtures/databases/atlas/atlas"
THREADS       = 16
merged_reads  = True
merged_prefix = "MergedReads-001"
metadata_path = None
bin_all       = True

def COLLECT_ALL_INPUT():
    INPUTS = []
    INPUTS.append("samples.tsv") #rule atlas_init
    if merged_reads:
        INPUTS.append("logs/concatReads.log")
    INPUTS.append("logs/formatSamples.log")
    INPUTS.append("finished_QC")
    return INPUTS

def COLLECT_INIT_INPUT():
    INPUTS = []
    INPUTS.append(fastq_dir)
    if merged_reads:
        INPUTS.append("logs/concatReads.log")
    return INPUTS

def COLLECT_FORMAT_ARGS():
    if bin_all:
        binarg = "--all"
    else:
        binarg = ""
    if metadata_path:
        return f"-s samples.tsv -m {metadata_path} {binarg}"
    else:
        return f"-s samples.tsv {binarg}"

# ---- Rule all: run all targets
rule all:
    input: COLLECT_ALL_INPUT()


# ---- Assembly Rules
rule concatReads:
    '''
    Use concatReads.py to merge read files  
    '''
    input: 
        fastq_dir = fastq_dir
    output: expand(os.path.join(fastq_dir, "{merged_prefix}_R{n}.fastq.gz"), merged_prefix = merged_prefix, n = [1,2])
    log: "logs/concatReads.log"
    benchmark: "benchmarks/concatReads.bmk"
    params:
        prefix = os.path.join(fastq_dir, f"{merged_prefix}"),
        d = 0,
        r = 0,
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
    shell:
        "scripts/concatReads.py -i {input} -o {params.prefix} -d {params.d} -r {params.r};"
        'echo Created merged reads files: {output} at {params.now} > {log}'


rule init_atlas:
    '''
    Initiate atlas by providing paths for
    Database directory and fastq files. 
    Outputs sample and config files.
    Should be run after concatReads. 
    ''' 
    input: COLLECT_INIT_INPUT()
    output:
        samples = "samples.tsv",
        config  = "config.yaml",
    log: "logs/init_atlas.log"
    benchmark: "benchmarks/init_atlas.bmk"
    params:
        fastq_dir    = fastq_dir,
        database_dir = database_dir
    threads: THREADS
    shell:
        "(atlas init --threads {threads} --db-dir {params.database_dir} {params.fastq_dir}) 2> {log}"


rule formatSamples:
    '''
    Use formatSamples.py to adjust samples.tsv file  
    '''
    input: "samples.tsv"
    output: "logs/formatSamples.log"
    log: "logs/formatSamples.log"
    benchmark: "benchmarks/formatSamples.bmk"
    params:
        args = COLLECT_FORMAT_ARGS(),
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
    shell:
        "python scripts/formatSamples.py {params.args};"
        'echo Modified bin groups in {input} at {params.now} > {log}'


rule atlas_qc:
    '''
    Run atlas qc. Provide path for
    Config file explicitly. 
    Outputs: finished_QC flag, 
    ref dir (bbmap), logs dir, reports dir, sample(s) dir(s), 
    Stats dir 
    ''' 
    input: 
        config = "config.yaml",
        format_proof = "logs/formatSamples.log"
    output: "finished_QC"
    log: "logs/atlas_qc.log"
    benchmark: "benchmarks/atlas_qc.bmk"
    threads: THREADS
    params: 
        mem = 128
    shell:
        "(atlas run qc -j {threads} --max-mem {params.mem} --config-file {input.config}) 2> {log};"
        "touch {output}"
