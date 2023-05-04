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

rule run_grist_gather:
    """
    Runs grist gather
    """
    input: "logs/create_grist_config_file.log"
    params:
        config = os.path.join(working_dir, "grist_config.yaml"),
        working_dir = working_dir
    output: os.path.join(working_dir, "finished_grist_gather")
    conda: 'grist'
    log: "logs/run_grist_gather.log"
    threads: int(config["threads"])
    shell:
        """
        cd {params.working_dir}
        (genome-grist run grist_config.yaml gather_reads -j {threads} ) 2> {log}
        touch {output}
        """

rule run_grist_repair:
    """
    Removes unavailable genomes from grist gather results
    """
    input: os.path.join(working_dir, "finished_grist_gather")
    params:
        working_dir = working_dir,
        script: os.path.join(config["parent_dir"], "workflows/scripts/repair_grist_gather_files.py")
    output: os.path.join(working_dir, "finished_grist_repair")
    conda: 'grist'
    log: "logs/run_grist_repair.log"
    shell:
        """
        cd {params.working_dir}
        cp {params.script} .
        python repair_grist_gather_files.py --grist_output_folder grist 2> {log}
        touch {output}
        """

rule run_grist:
    """
    Finishes the grist process
    """
    input: os.path.join(working_dir, "finished_grist_repair")
    params:
        config = os.path.join(working_dir, 'grist_config.yaml'),
        working_dir = working_dir
    output: os.path.join(working_dir, "finished_grist")
    conda: 'grist'
    log: "logs/run_grist.log"
    threads: int(config["threads"])
    shell:
        """
        cd {params.working_dir}
        (genome-grist run grist_config.yaml summarize_gather summarize_mapping -j {threads} -p) 2> {log}
        cd -
        touch {output}
        """

rule get_genomes_for_dereplication:
    """
    Get genomes from user or from grist
    """
    input:
        os.path.join(working_dir, "finished_grist"),
    params:
        user_genomes = config.get("genome")
        working_dir = working_dir,
    output: os.path.join(working_dir, "finished_genome_copying")
    run:
        shell("cd params.working_dir")
        shell("mkdir -p temp_genomes_folder")
        shell("cp grist/genomes/*fna.gz temp_genomes_folder")
        if params.user_genomes == "None" :
            continue
        else:
            shell("cp input.user_genomes/* temp_genomes_folder")
        shell("gunzip temp_genomes_folder/*gz")
        shell("touch output")


rule run_dereplicate_genomes:
    """
    Dereplicates genomes that come from Grist and from the user, by default there are no user-provided genomes
    """
    input: os.path.join(working_dir, "finished_genome_copying")
    params:
        working_dir = working_dir
    output: os.path.join(working_dir, "finished_genome_dereplication")# "finished dereplication"
    conda: 'drep'
    log: "logs/run_genome_replication.log"
    threads: int(config["threads"])
    shell:
    shell:
        """
        cd {params.working_dir}
        (dRep dereplicate derep_genomes/ -g temp_genomes_folder/ -p {threads}) 2> {log}
        cp derep_genomes/log/logger.log {log}
        touch {output}
        """


rule clean_after_grist:
    """
    Cleans folders after running grist, reorganized files
    """
    input: os.path.join(working_dir, "finished_genome_dereplication")
    params:
        working_dir = working_dir
    output: os.path.join(working_dir, "cleaned_after_grist")
    log: "logs/cleaned_grist.log"
    shell:
        """
        cd {params.working_dir}
        mkdir -p reference_genomes
        cp derep_genomes/dereplicated_genomes/* reference_genomes
        cd reference_genomes
        for i in `ls *fna|sed 's/.fna//g';do mv ${i}.fna ${i}.fasta;done
        cd -
        mv grist/reports logs
        mv grist/*trim .

        echo 'clean_after_grist' > {log}
        touch {output}
        """