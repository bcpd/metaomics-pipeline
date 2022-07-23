#!/usr/bin/python3

rule create_folder_structure_metagenomics:
    '''
    Creates and organized folder structure to store the important files
    after running Atlas QC, assembly, binning, and annotation
    ''' 
    input: 
        genomes_complete = os.path.join(working_dir, "finished_genomes"),
        refseeker_completed = os.path.join(working_dir, "finished_referenceseeker"),
        bakta_completed = os.path.join(working_dir, "finished_bakta"),
        dram_completed = os.path.join(working_dir, "finished_DRAM")
    output: os.path.join(working_dir, "logs/Creation_output_structure_metagenomics.log")
    log: os.path.join(working_dir, "logs/Creation_output_structure_metagenomics.log")
    params: 
        working_dir = working_dir
    shell:
        """
        cd {params.working_dir}
        mkdir metagenomics
        mkdir metagenomics/trimmed_reads
        mkdir metagenomics/assemblies
        mkdir metagenomics/MAGs
        mkdir metagenomics/taxonomic_annotations
        mkdir metagenomics/taxonomic_annotations/gtdb-tk
        mkdir metagenomics/taxonomic_annotations/referenceseeker
        mkdir metagenomics/MAGs/fasta
        mkdir metagenomics/MAGs/reports
        mkdir metagenomics/functional_annotations
        mkdir metagenomics/functional_annotations/GFF3
        """

rule reorganize_files_metagenomics:
    '''
    Copied important files from Atlas
    '''
    input:
        genomes_complete = os.path.join(working_dir, "finished_genomes"),
        refseeker_completed = os.path.join(working_dir, "finished_referenceseeker"),
        bakta_completed = os.path.join(working_dir, "finished_bakta"),
        dram_completed = os.path.join(working_dir, "finished_DRAM"),
        folders_created = os.path.join(working_dir, "logs/Creation_output_structure_metagenomics.log"),
    params:
        refseeker_file = os.path.join(working_dir, "refseeker.tsv"),
        bakta_file = os.path.join(working_dir, "bakta.tsv"),
    output: os.path.join(working_dir, "finished_metagenomics_cleanup")
    log: os.path.join(working_dir, "logs/Atlas_metagenomics_cleanup.log")
    benchmark: os.path.join(working_dir, "benchmarks/Atlas_metagenomics_cleanup.bmk")
    params:
        working_dir = working_dir
    shell:
        """
        cd {params.working_dir}
        # Copy quality-controlled reads
        cp */sequence_quality_control/*fastq.gz metagenomics/trimmed_reads/
        # Copy genome bins and related reports
        cp genomes/genomes/*fasta metagenomics/MAGs/fasta
        cp genomes/checkm/completeness.tsv metagenomics/MAGs/reports
        cp genomes/counts/* metagenomics/MAGs/reports

        # Copy assemblies, predicted genes, predicted proteins, and related annotations
        cp */assembly/*final_contigs.fasta metagenomics/assemblies/
        cp {params.bakta_file} metagenomics/functional_annotations
        cp genomes/annotations/genes/MAG*f?a metagenomics/functional_annotations
        cp DRAM/annotations/annotations.tsv  metagenomics/functional_annotations/dram_annotations.tsv
        cp DRAM/annotations/distill/product.tsv  metagenomics/functional_annotations/dram_product.tsv
        cp DRAM/annotations/distill/metabolism_summary.xlsx metagenomics/functional_annotations/dram_metabolism_summary.xlsx

        # Copying taxonomic annotations
        cp genomes/taxonomy/gtdb/classify/*summary.tsv metagenomics/taxonomic_annotations/gtdb-tk/
        cp {params.refseeker_file} metagenomics/taxonomic_annotations/referenceseeker

        # Copy stats from atlas process
        cp -r stats logs/
        cp -r logs/benchmarks/ benchmarks
        echo 'Copied metagenomics files to final folder' > {log}
        touch {output}
        """
