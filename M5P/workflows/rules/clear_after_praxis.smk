
import os
from pathlib import Path
import itertools
import pandas as pd


THREADS = config["threads"]
ALIGNER = "bt2"
METHOD = "salmon"
samples = "sample.tsv"
contrast = config["experimental_contrast"]
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


#rule deseq2:
#    input:
#        counts = os.path.join(working_dir, "salmon/counts.tsv")
#    output:
#        tables = os.path.join(working_dir, "salmon/DESeq2.tsv")
#    params:
#        samples = samples_table["SampleID"],
#        data = config["samples"],
#        contrasts = contrasts
#    conda:
#        "workflows/envs/deseq2.yaml"
#    script:
#        "workflows/scripts/deseq.R"


#rule report:
#    input:
#        counts = os.path.join(output_directory, "salmon/counts.tsv"),
#        tables = expand(DE_out, contrasts = contrasts)
#    output:
#        dereport_html = "DE_Report.html"
#    params:
#        data = config["samples"],
#        contrasts = contrasts,
#        ALIGNER = ALIGNER,
#        METHOD = METHOD,
#        PROJECT = project_name
#    script:
#        "workflows/scripts/DE_report.Rmd"