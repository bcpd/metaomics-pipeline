#!/usr/bin/python3

rule create_grist_config_file:
    """
    Creates a configurations file with the location of the grist database and the raw reads
    """
    input:
       samples_folder = fastq_dir2 # Metatranscriptomics reads
       #grist_database_dir = os.path.join(database_dir, "grist")
    params:
        grist_output_folder = os.path.join(working_dir, "grist"),
        max_memory = max_memory
    output: os.path.join(working_dir, "grist_config.yaml")
    conda:
       'M5P'
    log: os.path.join(working_dir, "log/create_grist_config_file.log")
    shell:
        "mkdir -p {params.grist_output_folder};"
        "python workflows/scripts/create_grist_config_file.py -s {input.samples_folder} -d ~/M5P_databases/grist -o {output} -m {max_memory} 2> {log}"

rule run_grist:
    '''
    Runs grist using raw sequences
    ''' 
    input: os.path.join(working_dir, "grist_config.yaml")
    output: os.path.join(working_dir, "grist/reports")
    conda: 'grist'
    log: os.path.join(working_dir, "logs/run_grist.log")
    threads: threads
    shell:
        "(genome-grist run {input} summarize_gather summarize_mapping -j {threads} -p) 2> {log}"