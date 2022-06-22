rule get_dram:
    """
    Installs DRAM in a docker image. This is to avoid issues with older versions of Ubuntu.    
    Once the program and the database are installed, it copies a log file to the databases directory.
    This log file (dram_setup_complete.log) is the output of the rule, so if it DRAM needs
    to be run again, we do not need to redownload the databases.
    """
    input:
        atlas_complete: os.path.join(working_dir, "finished_genomes")
    output:
        os.path.join(working_dir, "logs/dram_setup_complete.log")
    params:
        log_folder: os.path.join(working_dir, "logs/")
    log:
        os.path.join(working_dir, "logs/get_dram.log")
    shell:
        """
        docker pull continuumio/miniconda3
        docker run -i -d --name DRAM continuumio/miniconda3
        docker exec DRAM mkdir -p data out scripts logs 
        docker cp workflows/scripts/DRAM_setup.sh DRAM:/scripts
        docker exec DRAM /bin/bash /scripts/DRAM_setup.sh
        docker exec DRAM touch /logs/dram_setup_complete.log
        docker cp DRAM:/logs/dram_setup_complete.log {log_folder}
        docker cp DRAM:/logs/dram_setup_complete.log ~/M5P_databases
        echo 'Created dram container and setup databases' > {log}
        docker stop DRAM
        """

rule annotate_genomes:
    """
    Run DRAM inside the conda image, uses the predicted proteins as input
    """
    input:
       dram_setup_complete: os.path.join("~/M5P_databases/dram_setup_complete.log"),
       atlas_genome_complete: os.path.join(working_dir, "logs/Atlas_metagenomics_cleanup.log") 
    output: os.path.join(working_dir, "logs/DRAM_annotate.log")
    log: os.path.join(working_dir, "logs/DRAM_annotate.log")
    shell:
        """
        docker restart DRAM
        docker cp metagenomics/MAGs/fasta/*fasta DRAM:/genomes/
        docker cp workflows/scripts/annotate_genomes.sh DRAM:/scripts
        docker exec -t DRAM /bin/bash /scripts/annotate_genomes.sh 2> {log}
        """


rule copy_DRAM_annotations:
    """
    Copies the DRAM files from the docker image
    """
   input: os.path.join(working_dir, "logs/DRAM_annotate.log")
   output: os.path.join(working_dir, "logs/DRAM_copy_results.log") 
   log: os.path.join(working_dir, "logs/DRAM_copy_results.log")
   shell:
       """
       docker cp DRAM:out/ metagenomics/functional_annotations/"
       docker cp DRAM:logs/* /logs/"
       docker stop DRAM    
       touch {log}
       """
