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
        fasta = COLLECT_SALMON_REFERENCE()
    conda:
        'praxis'
    output:
        directory(output_directory + "reference/salmon_quasi")
    threads: THREADS
    run:
        shell("salmon index -t {input.fasta} -i {output} --type quasi -k 31")


rule salmon_quant:
    """
    Generate directories containing count files with salmon (quasi mode).
    """
    input:
        fastq_1= os.path.join(config["fastq_metatranscriptomics"], "{sample}_R1_001.fastq.gz"),
        fastq_2= os.path.join(config["fastq_metatranscriptomics"], "{sample}_R2_001.fastq.gz"),
        salmon_dir = output_directory + "reference/salmon_quasi"
    output:
        directory(os.path.join(output_directory, "/salmon/{sample}"))
    threads: THREADS
    conda:
        'praxis'
    log:
        "logs/salmon/{sample}.log"
    benchmark:
        "benchmarks/salmon/{sample}.bmk"
    run:
        shell("salmon quant -i {input.salmon_dir} -l A -p {threads} --validateMappings \
        -1 {input.fastq_1} -2 {input.fastq_2} -o {output} 2> {log}")

rule salmon_quant2:
    """
    Generate directories containing count files with salmon (quasi mode).
    """
    input:
        fastq_dir = fastq_metatranscriptomics,
        salmon_index = os.path.join(working_dir, "finished_salmon_index")
    output:
        directory(os.path.join(working_dir, "finished_salmon_quant"))
    params:
        working_dir = working_dir
    threads: threads
    shell:
        """
        for i in {input.folder}/*fastq.gz;do cp $i {working_dir};done
        cd {params.working_dir}
        for i in `ls *R1_001.fastq.gz|sed 's/_R1_001.fastq.gz//g'`;do \
        salmon quant -i reference/salmon_quasi -l A -p 10 --validateMappings -1 ${{i}}_R1_001.fastq.gz -2 ${{i}}_R2_001.fastq.gz -o salmon/${{i}} ;done
        """


rule salmon_quant_table:
    """
    Generate a count table with salmon.
    """
    input:
        expand(os.path.join(output_directory, "salmon/{sample}"), sample = sample)
    conda:
        'praxis'
    output:
        os.path.join(output_directory, "salmon/counts.tsv")
    run:
        shell("salmon quantmerge --quants {input} --column numreads -o {output}")


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
