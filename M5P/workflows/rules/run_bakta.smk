rule run_bakta:
    '''
    Run bakta annotation.
    Requires database directory, output directory, and input directory containing MAGs.
    Outputs one annotation file for all the MAGs.
    '''
    input:
        binning_completed = os.path.join(working_dir, "finished_binning"),
        genomes_folder = os.path.join(working_dir, "genomes/genomes")
    output: os.path.join(working_dir, "logs/run_bakta.log")
    benchmark: os.path.join(working_dir, "benchmarks/bakta.bmk")
    conda:
        'bakta'
    log: os.path.join(working_dir, "logs/run_bakta.log")
    params:
        working_dir  = working_dir
    conda:
        'bakta'
    threads: config["threads"]
    shell:
        """
        cd {params.working_dir}
        mkdir -p bakta
        cd bakta
        cp -r ../genomes/genomes/ .
        cd genomes
        for i in *.fasta; do bakta --db ~/M5P_databases/bakta/db $i --output .. --threads {threads}; done
        cd ..
        echo SRA$'\t'Sequence_Id$'\t'Type$'\t'start_position$'\t'end_position$'\t'Strand$'\t'Locus_Tag$'\t'Gene$'\t'Product$'\t'Product$'\t'DbXrefs > bakta.tsv
        for f in *.tsv; do
        if [[ $f != *"hypotheticals"* ]] ; then
            i=$(basename $f .tsv)
            sed -e '1,3d' $f | sed "s/^/\t$i /" > $i.filename.tsv
            cat *filename* >> bakta.tsv
        fi
        done
        sed -i 's/^\t//g' bakta.tsv
        rm *filename*
        mv bakta.tsv {params.working_dir}
        cd {params.working_dir}
        echo 'bakta completed' > {log}
        """