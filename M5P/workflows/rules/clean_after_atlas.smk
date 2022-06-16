#!/usr/bin/python3

rule create_folder_structure_metagenomics:
    '''
    Creates and organized folder structure to store the important files
    after running Atlas QC, assembly, binning, and annotation
    ''' 
    input: 
        atlas_genecatalog_log = os.path.join(working_dir, "logs/atlas_genecatalog.log")
        atlas_gtdbtk = os.path.join(working_dir, "genomes/taxonomy/gtdb.log")
        refseeker_file = os.path.join(working_dir, "referenceseeker.tsv)
    output: os.path.join(working_dir, "logs/Creation_output_structure_metagenomics.log")
    log: os.path.join(working_dir, "logs/Creation_output_structure_metagenomics.log")
    params: 
        working_dir = working_dir,
    shell:
        """
        cd {working_dir}
        mkdir metagenomics metagenomics/trimmed_reads
        mkdir metagenomics/MAGs metagenomics/MAGS/reports metagenomics/MAGS/fasta
        mkdir metagenomics/functional_annotations metagenomics/functional_annotations/GFF3
        mkdir metagenomics/taxonomic_annotations
        mkdir metagenomics/taxonomic_annotations/referenceseeker
        mkdir metagenomics/taxonomic_annotations/gtdb-tk/
        mkdir metagenomics/assemblies
        """

rule create_folder_structure_metatranscriptomics:
    '''
    Creates and organized folder structure to store the important files
    after running praxis
    ''' 
    input: 
        config = os.path.join(working_dir, "logs/praxis.log")
    output: os.path.join(working_dir, "logs/Creation_output_structure_metatranscriptomics.log")
    params: 
        working_dir = working_dir,
    shell:
        """
        cd {working_dir}
        mkdir metatranscriptomics
        mkdir metatranscriptomics/trimmed_reads
        mkdir metatranscriptomics/grist
        mkdir metatranscriptomics/transcript_counts
        mkdir metatranscriptomics/interleaved_reads
        """


rule reorganize_files_metagenomics:
    '''
    Copied important files from Atlas 
    ''' 
    input: 
        gene_catalog_log = os.path.join(working_dir, "logs/Creation_output_structure_metagenomics.log")
        atlas_gtdbtk_log = os.path.join(working_dir, "genomes/taxonomy/gtdb.log")
        refseeker_file = os.path.join(working_dir, "refseeker.tsv")
    output: os.path.join(working_dir, "logs/Atlas_metagenomics_cleanup.log")
    log: os.path.join(working_dir, "logs/Atlas_metagenomics_cleanup.log")
    benchmark: os.path.join(working_dir, "benchmarks/Atlas_metagenomics_cleanup.bmk")
    params: 
        working_dir = working_dir
    shell:
        """
        cd {working_dir}
        # Copy quality-controlled reads	
        cp */sequence_quality_control/*fastq.gz metagenomics/trimmed_reads/
        cp genome/checkm/* metagenomics/MAGs/fasta
        # Copy genome bins and related reports
        cp */binning/DASTool/*.eval  metagenomics/MAGs/reports
        cp */binning/DASTool/*DASTool_summary.tsv metagenomics/MAGs/reports
        # Copy assemblies, predicted genes, predicted proteins, and related annotations
        cp */assembly/*final_contigs.fasta metagenomics/assemblies/
        cp */annotation/predicted_genes/*gff metagenomics/functional_annotations/GFF3
        cp Genecatalog/*f?a metagenomics/functional_annotations
        cp Genecatalog/counts/ metagenomics/functional_annotations
        gunzip metagenomics/functional_annotations/*gz 
        # Copying taxonomic annotations
        cp genomes/taxonomy/gtdb/classify/*summary.tsv metagenomics/taxonomic_annotations/gtdb-tk/
        cp {refseeker_file} metagenomics/taxonomic_annotations/referenceseeker
        # Copy stats from atlas process
        cp -r stats logs/
        cp logs/benchmarks/*bmk benchmarks
        echo 'Copied metagenomics files to final folder' > {log}
        """


rule reorganize_files_transcriptomic:
    '''
    Copied important files from Praxis process 
    ''' 
    input: os.path.join(working_dir, "logs/Creation_output_structure_metatranscriptomics.log")
    output: os.path.join(working_dir, "logs/Praxis_cleanup.log")
    log: os.path.join(working_dir, "logs/Praxis_cleanup.log")
    benchmark: os.path.join(working_dir, "benchmarks/Praxis_cleanup.bmk")
    params: 
        working_dir = working_dir,
    shell:
        """
        cd {working_dir}
        # Copy files
        echo 'Copied metatranscriptomics files to final folder' > {log}
        """