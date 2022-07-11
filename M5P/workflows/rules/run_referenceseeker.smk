rule run_referenceseeker:
    '''
    Run referenceseeker for taxonomic annotation.
    Requires database directory, output directory, and input directory containing MAGs.
    Outputs one annotation file for all the MAGs and individually annotated MAGs.
    '''
    input: os.path.join(working_dir, "finished_binning")
    output: os.path.join(working_dir, "metagenomics/taxonomic_annotations/refseeker.tsv")
    benchmark: os.path.join(working_dir, "benchmarks/refseeker.bmk")
    conda:
        'referenceseeker'
    log: os.path.join(working_dir, "log/run_referenceseeker.log")
    params:
        working_dir  = working_dir
    conda:
        'referenceseeker'
    threads: THREADS
    shell:
        "(workflows/scripts/run_referenceseeker.sh -i {input} -d ~/M5P_databases/referenceseeker -o {params.working_dir}) 2> {log}"
