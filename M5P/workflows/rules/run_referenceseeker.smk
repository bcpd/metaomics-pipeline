rule run_referenceseeker:
    '''
    Run referenceseeker for taxonomic annotation.
    Requires database directory, output directory, and input directory containing MAGs.
    Outputs one annotation file for all the MAGs and individually annotated MAGs.
    '''
    input: os.path.join(working_dir, "finished_binning")
    output: os.path.join(working_dir, "refseeker.tsv")
    benchmark: os.path.join(working_dir, "benchmarks/refseeker.bmk")
    log: os.path.join(working_dir, "log/run_refseeker.log")
    params:
        output_dir    = output_dir,
        #database_dir = database_dir
        working_dir  = working_dir
    threads: THREADS
    shell:
        "(workflows/scripts/run_referenceseeker.sh -i {input} -d ~/M5P_databases/referenceseeker -o {params.output_dir}) 2> {log}"
