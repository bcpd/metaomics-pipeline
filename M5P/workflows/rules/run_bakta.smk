rule run_bakta:
    '''
    Run bakta annotation.
    Requires database directory, output directory, and input directory containing MAGs.
    Outputs one annotation file for all the MAGs.
    '''
    input:
        binning_completed = os.path.join(working_dir, "finished_binning"),
        genomes_folder = os.path.join(working_dir, "genomes/genomes")
    output: os.path.join(working_dir, "bakta.tsv")
    benchmark: os.path.join(working_dir, "benchmarks/bakta.bmk")
    conda:
        'bakta'
    log: os.path.join(working_dir, "logs/run_bakta.log")
    params:
        working_dir  = working_dir,
        script = os.path.join(config["parent_dir"], "workflows/scripts/run_bakta.sh")
    conda:
        'bakta'
    threads: THREADS
    shell:
        "({params.script} -i {input.genomes/genomes} -d ~/M5P_databases/bakta -o {params.working_dir}) 2> {log}"
