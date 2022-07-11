import os
from pathlib import Path
import itertools
import pandas as pd


THREADS = config["THREADS"]
ALIGNER = config["ALIGNER"]
METHOD = config["METHOD"]

rule salmon_index:
    """
    Index a reference with salmon.
    """
    input:
        fasta = COLLECT_SALMON_REFERENCE()
    output:
        directory(output_directory + "reference/salmon_quasi"),
    threads: THREADS
    run:
        shell("salmon index -t {input.fasta} -i {output} --type quasi -k 31")


rule salmon_quant:
    """
    Generate directories containing count files with salmon (quasi mode).
    """
    input:
        fastq_1= os.path.join(config["fastq_dir2"], "{sample}_R1.fastq.gz",
        fastq_2= os.path.join(config["fastq_dir2"], "{sample}_R2.fastq.gz",
        salmon_dir = output_directory + "reference/salmon_quasi"
    output:
        directory(os.path.join(working_dir, "/salmon/{sample}"))
    threads: THREADS
    log:
        "log/salmon/{sample}.log"
    benchmark:
        "benchmarks/salmon/{sample}.bmk"
    run:
        shell("salmon quant -i {input.salmon_dir} -l A -p {threads} --validateMappings \
        -1 {input.fastq_1} -2 {input.fastq_2} -o {output} 2> {log}")


rule salmon_quant_table:
    """
    Generate a count table with salmon.
    """
    input:
        expand(os.path.join(working_dir, "/salmon/{sample}"), sample = sample)
    output:
        os.path.join(working_dir, "/salmon/counts.tsv")
    run:
        shell("salmon quantmerge --quants {input} --column numreads -o {output}")


rule deseq2:
    input:
        counts = os.path.join(working_dir, "/salmon/counts.tsv")
    output:
        tables = os.path.join(working_dir, "/salmon/DESeq2.tsv")
    params:
        samples = samples_table["SampleID"],
        data = config["samples"],
        contrasts = contrasts
    conda:
        "/envs/deseq2.yaml"
    script:
        "../scripts/deseq.R"


rule report:
    input:
        counts = os.path.join(working_dir, "/salmon/counts.tsv")
        tables = expand(DE_out, contrasts = contrasts),
    output:
        dereport_html = "DE_Report.html"
    params:
        data = config["samples"],
        contrasts = contrasts,
        ALIGNER = config["ALIGNER"],
        METHOD = config["METHOD"],
        PROJECT = project_name,
    script:
        "../scripts/DE_report.Rmd"
