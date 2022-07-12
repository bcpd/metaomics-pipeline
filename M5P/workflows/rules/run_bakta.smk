rule run_bakta:
    '''
    Run bakta annotation.
    Requires database directory, output directory, and input directory containing MAGs.
    Outputs one annotation file for all the MAGs.
    '''
    input: os.path.join(working_dir, "finished_binning")
    output: os.path.join(working_dir, "bakta.tsv")
    benchmark: os.path.join(working_dir, "benchmarks/bakta.bmk")
    conda:
        'bakta'
    log: os.path.join(working_dir, "log/run_bakta.log")
    params:
        working_dir  = working_dir,
        script = os.path.join(config["parent_dir"], "workflows/scripts/run_bakta.sh")
    conda:
        'bakta'
    threads: THREADS
    shell:
        "(/bin/bash {params.script} -i {input} -d ~/M5P_databases/bakta -o {params.working_dir}) 2> {log}"
