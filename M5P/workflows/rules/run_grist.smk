#!/usr/bin/python3

rule create_grist_config_file:
    """
    Creates a configurations file with the location of the grist database and the raw reads
    """
    input:
       samples_folder = fastq_metatranscriptomics # Metatranscriptomics reads
    params:
        grist_output_folder = os.path.join(working_dir, "grist"),
        max_memory = config["max_memory"],
        script = os.path.join(config["parent_dir"], "workflows/scripts/create_grist_config_file.py")
    output: os.path.join(working_dir, "grist_config.yaml")
    conda: 'M5P_test'
    log: os.path.join(working_dir, "log/create_grist_config_file.log")
    shell:
        "mkdir -p {params.grist_output_folder};"
        "python {params.script} -s {input.samples_folder} -d ~/M5P_databases/grist -o {output} -m {max_memory} 2> {log}"

rule run_grist:
    '''
    Runs grist using raw sequences
    ''' 
    input: os.path.join(working_dir, "grist_config.yaml")
    output: directory(os.path.join(working_dir, "grist/reports"))
    conda: 'grist'
    log: os.path.join(working_dir, "logs/run_grist.log")
    threads: THREADS
    shell:
        "(genome-grist run {input} summarize_gather summarize_mapping -j {threads} -p) 2> {log}"
