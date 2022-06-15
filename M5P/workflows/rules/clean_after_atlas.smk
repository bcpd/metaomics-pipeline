rule create_folder_structure:
    '''
    Creates and organized folder structure to store the important files
    after running Atlas QC, assembly, binning, and annotation
    ''' 
    input: 
        config = os.path.join(working_dir, "M5P_config.yaml"),
    log: os.path.join(working_dir, "logs/Atlas_cleanup.log")
    benchmark: os.path.join(working_dir, "benchmarks/Atlas_cleanup.bmk")
    params: 
        working_dir = working_dir,
        mem = 128
    shell:
        """
        cd working_dir
        mkdir metagenomics metagenomics/trimmed_reads
        mkdir metagenomics/MAGs metagenomics/MAGS/reports
        mkdir metagenomics/functional_annotations
        mkdir metagenomics/taxonomic_annotations
        mkdir metatranscriptomics metatranscriptomics/trimmed_reads
        mkdir metatranscriptomics/grist
        metatranscriptomics/transcript_counts
        metatranscriptomics/interleaved_reads
        mkdir logs benchmarks
        """

rule reorganize_files_metagenomics:
    '''
    Create organized folder structure to store the important files
    after running Atlas QC, assemlby. 
    ''' 
    input: 
        config = os.path.join(working_dir, "M5P_config.yaml"),
    log: os.path.join(working_dir, "logs/Atlas_cleanup.log")
    benchmark: os.path.join(working_dir, "benchmarks/Atlas_cleanup.bmk")
    params: 
        working_dir = working_dir,
        mem = 128
    shell:
        """
        cd {working_dir}
        mkdir metagenomics metagenomics/trimmed_reads
        mkdir metagenomics/MAGs metagenomics/MAGS/reports
        mkdir metagenomics/functional_annotations
        mkdir metagenomics/taxonomic_annotations
        mkdir metatranscriptomics metatranscriptomics/trimmed_reads
        mkdir metatranscriptomics/grist
        metatranscriptomics/transcript_counts
        metatranscriptomics/interleaved_reads
        mkdir logs benchmarks
        """