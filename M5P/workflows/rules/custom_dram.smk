#!/usr/bin/python3

rule get_dram:
    '''
    Installs DRAM in a docker image. This is to avoid issues with older versions of Ubuntu.    
    Once the program and the database are installed, it copies a log file to the databases directory.
    This log file (dram_setup_complete.log) is the output of the rule, so if it DRAM needs
    to be run again, we do not need to redownload the databases.
    '''
    input:
        atlas_complete= os.path.join(working_dir, "finished_genomes")
    output:
        os.path.join(working_dir, "logs/dram_setup_complete.log")
    params:
        log_folder = os.path.join(working_dir, "logs/"),
        script = os.path.join(config["parent_dir"], "workflows/scripts/DRAM_setup.sh")
    log:
        os.path.join(working_dir, "logs/get_dram.log")
    shell:
        """
        docker pull continuumio/miniconda3
        docker run -i -d --name DRAM continuumio/miniconda3
        docker exec DRAM mkdir -p data out scripts logs genomes proteins
        docker cp {params.script} DRAM:/scripts
        docker exec DRAM /bin/bash /scripts/DRAM_setup.sh
        docker exec DRAM touch /logs/dram_setup_complete.log
        docker cp DRAM:/logs/dram_setup_complete.log {log_folder}
        docker cp DRAM:/logs/dram_setup_complete.log ~/M5P_databases
        echo 'Created DRAM container and setup databases' > {log}
        docker stop DRAM
        """

rule annotate_genomes:
    '''
    Run DRAM inside the conda image, uses the predicted proteins as input
    '''
    input:
       dram_setup_complete = os.path.join("~/M5P_databases/dram_setup_complete.log"),
       atlas_genome_complete = os.path.join(working_dir, "logs/atlas_genomes.log")
    output: os.path.join(working_dir, "logs/DRAM_annotate.log")
    log: os.path.join(working_dir, "logs/DRAM_annotate.log")
    shell:
        """
        docker restart DRAM
        docker exec DRAM rm /genomes/*
        for i in genome/annotations/genes/MAG*faa; do docker cp $i DRAM:/genomes;done
        docker cp workflows/scripts/DRAM_annotate_proteins.sh DRAM:/scripts
        docker exec -t DRAM /bin/bash /scripts/DRAM_annotate_proteins.sh 2> {log}
        """


rule copy_DRAM_annotations:
    '''
    Copies the DRAM output files from the docker image to the local directory
    '''
    input: os.path.join(working_dir, "logs/DRAM_annotate.log")
    output: os.path.join(working_dir, "logs/DRAM_copy_results.log")
    params:
        output_folder= os.path.join(working_dir, "DRAM")
    log: os.path.join(working_dir, "logs/DRAM_copy_results.log")
    shell:
        """
        mkdir {params.output_folder}
        docker cp DRAM:out/annotations {params.output_folder}
        docker cp DRAM:logs/* /logs/
        docker stop DRAM    
        touch {log}
        """
