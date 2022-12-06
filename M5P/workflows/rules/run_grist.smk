#!/usr/bin/python3

rule prepare_grist_folders:
    """
    Creates folders in a strucutre that grist can use
    """
    input:
        samples_folder  = fastq_metatranscriptomics
    params:
        grist_output_folder = os.path.join(working_dir, "grist"),
        working_dir = working_dir,
        temp_folder = os.path.join(working_dir, "temp"),
    output: os.path.join(working_dir, "logs/create_grist_folders.log")
    log: "logs/create_grist_folders.log"
    shell:
        """
        mkdir -p {params.grist_output_folder}
        mkdir -p {params.temp_folder}
        cp -r {input.samples_folder} {params.temp_folder}
        cd {params.working_dir}
        for i in temp/*;do mv $i input;done
        touch {log}
        """


rule create_grist_config_file:
    """
    Creates a configurations file with the location of the grist database and the raw reads
    """
    input:
        os.path.join(working_dir, "logs/create_grist_folders.log")
    params:
        working_dir = working_dir,
        max_memory = config["max_memory"],
        script = os.path.join(config["parent_dir"], "workflows/scripts/create_grist_config_file.py"),
        config_file = os.path.join(working_dir, "grist_config.yaml")
    output: os.path.join(working_dir, "logs/create_grist_config_file.log")
    conda: 'M5P'
    log: "logs/create_grist_config_file.log"
    shell:
        """
        cp {params.script} .
        python create_grist_config_file.py -d ~/M5P_databases/grist -s input -o grist_config.yaml -m {params.max_memory}
        touch {log}
        """

rule run_grist:
    '''
    Runs grist using raw sequences
    '''
    #input: os.path.join(working_dir, "grist_config_created")
    input: "logs/create_grist_config_file.log"
    params:
        config = os.path.join(working_dir, 'grist_config.yaml'),
        working_dir = working_dir
    output: os.path.join(working_dir, "finished_grist")# "finished_grist"
    conda: 'grist'
    log: "logs/run_grist.log"
    threads: config["threads"]
    shell:
        """
        cd {params.working_dir}
        (genome-grist run grist_config.yaml summarize_gather summarize_mapping -j {threads} -p) 2> {log}
        cd -
        touch {output}
        """


rule clean_after_grist:

    """
    Cleans folders after running grist, reorganized files
    """
    input: os.path.join(working_dir, "finished_grist")
    params:
        working_dir = working_dir
    output: os.path.join(working_dir, "cleaned_after_grist")
    log: "logs/cleaned_grist.log"
    shell:
        """
        cd {params.working_dir}
        mv grist/genomes .
        gunzip genomes/*gz
        mv genomes reference_genomes
        mv grist/reports logs
        mv grist/abundtrim .
        rm -fr grist temp create_grist_config_file.py input  genbank_cache
        echo 'clean_after_grist' > {log}
        cd -
        touch {output}
        """