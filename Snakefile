#!/usr/bin/python3


# Module Imports
#----------------------------------------------------------------------------#

include: "rules/common.smk"

import re
import os
import sys
from os import listdir
from os.path import isfile, join
import itertools
import pandas as pd
from pathlib import Path
import argparse
import psutil
import copy
from datetime import datetime
now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
import yaml
import configparser


# Configuration 
#----------------------------------------------------------------------------#

EXPERIMENT_TYPE = config["EXPERIMENT_TYPE"]

MAG_ANNOTATION = config["MAG_ANNOTATION"]

ASSEMBLY_STRATEGY = config["ASSEMBLY_STRATEGY"]

GENOME_REFERENCES = config["GENOME_REFERENCES"]

ANNOTATION_DATABASES = config["ANNOTATION_DATABASES"]

# Experimental design


samples_table = pd.read_csv(config["samples"], sep="\t")
sampleID_list = samples_table["SampleID"]


EXP_METADATA = config["METADATA"]
EXP_DESIGN = config["EXP_DESIGN"]
EXP_CONTRAST = config["EXP_CONTRAST"]

# Computational resources
N_CORES = config["N_CORES"]
MAX_MEMORY = config["MAX_MEMORY"]

#

METAGENOMICS_INPUT_FOLDER = config["MAX_MEMORY"]
METATRANSCRIPTOMICS_INPUT_FOLDER = config["MAX_MEMORY"]

OUTPUT_FOLDER = config["OUTPUT_FOLDER"]



# Begin Snakemake 


# ---- Targets rules
include: "workflow/rules/targets/targets.smk"

# ---- Quality control rules
include: "workflow/rules/qc/atlas_init.smk"
include: "workflow/rules/qc/atlas_qc.smk"

# ---- Assembly rules
include: "workflow/rules/assembly/atlas_assembly.smk"
include: "workflow/rules/assembly/atlas_binning.smk"


# ---- Contig annotation rules
include: "workflow/rules/annotation/atlas_annotation.smk"
include: "workflow/rules/annotation/dram_annotation.smk"
include: "workflow/rules/annotation/bakta_annotations.smk"
#include: "workflow/rules/annotation/integrate_annotations.smk"


# ---- Metatranscriptomics rules
  
include: "workflow/rules/metatranscriptomics/praxis_quality_control.smk"

 # Use this option is a link to an NCBI genome is provided
if GENOME_REFERENCES == "referemnce_url":
    include: "workflow/rules/metatranscriptomics/praxis_download.smk"
    annotate: false
elif GENOME_REFERENCES == "referemnce_file":     
      annotate: false
# Use this option if you want to make a denove metranscriptome assembly
elif GENOME_REFERENCES == "assemble_metatranscriptome":  
    include: "workflow/rules/metatranscriptomics/praxis_assemble.smk"
# Use this option if you want to find close relatives with grist
else GENOME_REFERENCES == "find_relatives":
    include: "workflow/rules/metatranscriptomics/grist.smk"
     annotate: false

include: "workflow/rules/metatranscriptomics/praxis_index_align.smk"
include: "workflow/rules/metatranscriptomics/praxis_quantify.smk"
include: "workflow/rules/metatranscriptomics/praxis_transcript_annotate.smk"


# ---- Rule all: run all targets
rule all:
    input: TARGET_ALL

