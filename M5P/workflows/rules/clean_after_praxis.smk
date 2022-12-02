#!/usr/bin/python3

rule create_folder_structure_metatranscriptomics:
    '''
    Creates and organized folder structure to store the important files
    after running praxis
    ''' 
    input: os.path.join(working_dir, "logs/salmon_quant_table.log")
    output: os.path.join(working_dir, "logs/Creation_output_structure_metatranscriptomics.log")
    params: 
        working_dir = working_dir
    log: os.path.join(working_dir, "logs/Creation_output_structure_metatranscriptomics.log")
    shell:
        """
        cd {params.working_dir}
        mkdir -p metatranscriptomics
        mkdir -p metatranscriptomics/trimmed_reads
        mkdir -p metatranscriptomics/reference_genomes
        mkdir -p metatranscriptomics/transcript_counts
        mkdir -p metatranscriptomics/logs
        touch {log}
        """


rule reorganize_files_transcriptomic:
    '''
    Copied important files from Praxis process
    ''' 
    input:
        creation_log = os.path.join(working_dir, "logs/Creation_output_structure_metatranscriptomics.log"),
        praxis_log = (working_dir, "logs/praxis.log")
    output: os.path.join(working_dir, "logs/Praxis_cleanup.log")
    log: os.path.join(working_dir, "logs/Praxis_cleanup.log")
    benchmark: os.path.join(working_dir, "benchmarks/Praxis_cleanup.bmk")
    params: 
        working_dir = working_dir
    shell:
        """
        cd {params.working_dir}
        # Copy files
        mv abundtrim/* metatranscriptomics/trimmed_reads
        cp -r reference_genomes/* metatranscriptomics/reference_genomes
        cp -r salmon/* metatranscriptomics/transcript_counts
        echo 'Copied metatranscriptomics files to final folder' > logs/Praxis_cleanup.log
        cp -r  logs/* metatranscriptomics/logs
        """