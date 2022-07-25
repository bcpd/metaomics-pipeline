#!/usr/bin/python3

import os
from pathlib import Path
import itertools
import pandas as pd


THREADS = config["threads"]
ALIGNER = "bt2"
METHOD = "salmon"
#contrasts = config["experimental_contrast"]
#design = config["experimental_design"]
output_directory = config["working_dir"]
project_name = "Microbial mixtures"

rule salmon_index:
    """
    Index a reference with salmon.
    """
    input:
        dram_output = os.path.join(working_dir, "finished_DRAM_annotate_reference_genomes")
    params:
        index = os.path.join(working_dir, "reference/salmon_quasi"),
        fasta = os.path.join(working_dir, "reference_genomes/annotations/genes.fna")
    output: os.path.join(working_dir, "finished_salmon_index")
    threads: THREADS
    shell:
        "salmon index -t {params.fasta} -i {params.index} --type quasi -k 31; touch {output}"


rule salmon_quant2:
    """
    Generate directories containing count files with salmon (quasi mode).
    """
    input:
        fastq_dir = config["fastq_metatranscriptomics"],
        salmon_index = os.path.join(working_dir, "finished_salmon_index")
    output:
        os.path.join(working_dir, "finished_salmon_quant")
    params:
        working_dir = working_dir
    threads: THREADS
    shell:
        """
        for i in {input.fastq_dir}/*fastq.gz;do cp $i {params.working_dir};done
        cd {params.working_dir}
        for i in `ls *R1_001.fastq.gz|sed 's/_R1_001.fastq.gz//g'`;do \
        salmon quant -i reference/salmon_quasi -l A -p {threads} --validateMappings -1 ${{i}}_R1_001.fastq.gz -2 ${{i}}_R2_001.fastq.gz -o salmon/${{i}} ;done
        touch finished_salmon_quant
        """


rule salmon_quant_table:
    """
    Generate a count table with salmon.
    """
    input:
        os.path.join(working_dir, "finished_salmon_quant")
    params:
        working_dir = working_dir
    log:
        "logs/salmon_quant_table.log"
    output:
        os.path.join(working_dir, "salmon/counts.tsv")
    shell:
        """
        cd {params.working_dir}
        salmon quantmerge --quants salmon/* --column numreads -o salmon/counts.tsv > {log}
        """


rule deseq2:
    input:
        counts = os.path.join(output_directory, "salmon/counts.tsv")
    output:
        tables = os.path.join(output_directory, "salmon/DESeq2.tsv")
    params:
        samples = samples_table["SampleName"],
        data = config["samples"],
        #contrasts = contrasts
        #design = design
    conda:
        "workflows/envs/deseq2.yaml"
    script:
        "workflows/scripts/deseq.R"


rule report:
    """
    Required R libraries
    pander, ggplot2, stringr, EnhancedVolcano, reshape2, vegan
    """
    input:
        counts = os.path.join(output_directory, "salmon/counts.tsv"),
        tables = expand(DE_out, contrasts = contrasts)
    output:
        dereport_html = "DE_Report.html"
    conda:
        'praxis'
    params:
        data = config["samples"],
        #contrasts = contrasts,
        ALIGNER = ALIGNER,
        METHOD = METHOD,
        PROJECT = project_name
    script:
        "workflows/scripts/DE_report.Rmd"
