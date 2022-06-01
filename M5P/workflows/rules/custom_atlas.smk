#!/usr/bin/python3

#----------------------------------------------------------------------------#
#
# Microbiome Insights | Snakemake pipeline for microbial characterization
#
#----------------------------------------------------------------------------#


# ---- Atlas Assembly Rules
rule concatReads:
    '''
    Use concatReads.py to merge read files  
    '''
    input: 
        fastq_dir = fastq_dir
    output: expand(os.path.join(fastq_dir, "{merged_prefix}_R{n}.fastq.gz"), merged_prefix = merged_prefix, n = [1,2])
    log: os.path.join(working_dir, "logs/concatReads.log")
    benchmark: os.path.join(working_dir, "benchmarks/concatReads.bmk")
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
        samples = os.path.join(working_dir, "samples.tsv"),
        config  = os.path.join(working_dir, "config.yaml"),
    log: os.path.join(working_dir, "logs/init_atlas.log")
    benchmark: os.path.join(working_dir, "benchmarks/init_atlas.bmk")
    params:
        fastq_dir    = fastq_dir,
        database_dir = database_dir,
        working_dir  = working_dir
    threads: THREADS
    shell:
        "(atlas init --threads {threads} -w {params.working_dir} --db-dir {params.database_dir} {params.fastq_dir}) 2> {log}"


rule formatSamples:
    '''
    Use formatSamples.py to adjust samples.tsv file  
    '''
    input: os.path.join(working_dir, "samples.tsv")
    output: os.path.join(working_dir, "logs/formatSamples.log")
    log: os.path.join(working_dir, "logs/formatSamples.log")
    benchmark: os.path.join(working_dir, "benchmarks/formatSamples.bmk")
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
        config = os.path.join(working_dir, "config.yaml"),
        format_proof = os.path.join(working_dir, "logs/formatSamples.log")
    output: os.path.join(working_dir, "finished_QC")
    log: os.path.join(working_dir, "logs/atlas_qc.log")
    benchmark: os.path.join(working_dir, "benchmarks/atlas_qc.bmk")
    threads: THREADS
    params: 
        working_dir = working_dir,
        mem = 128
    shell:
        "(atlas run qc -j {threads} -w {params.working_dir} --max-mem {params.mem} --config-file {input.config}) 2> {log};"
        "touch {output}"


rule atlas_assembly:
    '''
    Run atlas assembly. Provide path for
    Config file explicitly. 
    Ensure QC step is complete.
    Creates <sample name>_contigs.fasta files in
    Sample dirs 
    Outputs: finished_assembly
    ''' 
    input: 
        config = os.path.join(working_dir, "config.yaml"),
        qc_proof = os.path.join(working_dir, "finished_QC")
    output: os.path.join(working_dir, "finished_assembly")
    log: os.path.join(working_dir, "logs/atlas_assembly.log")
    benchmark: os.path.join(working_dir, "benchmarks/atlas_assembly.bmk")
    threads: THREADS
    params: 
        working_dir = working_dir,
        mem = 128
    shell:
        "(atlas run assembly -j {threads} -w {params.working_dir} --max-mem {params.mem} --config-file {input.config}) 2> {log};"
        "touch {output}"



rule atlas_binning:
    '''
    Run atlas binning. Provide path for
    Config file explicitly. 
    Ensure assembly step is complete.
    Sample dirs 
    Outputs: finished_bining
    ''' 
    input: 
        config = os.path.join(working_dir, "config.yaml"),
        assembly_proof = os.path.join(working_dir, "finished_assembly")
    output: os.path.join(working_dir, "finished_binning")
    log: os.path.join(working_dir, "logs/atlas_binning.log")
    benchmark: os.path.join(working_dir, "benchmarks/atlas_binning.bmk")
    threads: THREADS
    params: 
        working_dir = working_dir,
        mem = 128
    shell:
        "(atlas run binning -j {threads} -w {params.working_dir} --max-mem {params.mem} --config-file {input.config}) 2> {log};"
        "touch {output}"


rule atlas_genecatalog:
    '''
    Run atlas annotation. Provide path for
    Config file explicitly. 
    Ensure binning step is complete.
    Outputs: annotation files
    ''' 
    input: 
        config = os.path.join(working_dir, "config.yaml"),
        assembly_proof = os.path.join(working_dir, "finished_binning")
    output: os.path.join(working_dir, "finished_genecatalog")
    log: os.path.join(working_dir, "logs/atlas_genecatalog.log")
    benchmark: os.path.join(working_dir, "benchmarks/atlas_genecatalog.bmk")
    threads: THREADS
    params: 
        working_dir = working_dir,
        mem = 128
    shell:
        "(atlas run genecatalog -j {threads} -w {params.working_dir} --max-mem {params.mem} --config-file {input.config}) 2> {log};"
        "touch {output}"
