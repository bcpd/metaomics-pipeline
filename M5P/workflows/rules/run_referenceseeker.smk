#!/usr/bin/python3

rule run_referenceseeker:
    '''
    Run referenceseeker for taxonomic annotation.
    Requires database directory, output directory, and input directory containing MAGs.
    Outputs one annotation file for all the MAGs and individually annotated MAGs.
    '''
    input:
        genomes_folder = os.path.join(working_dir, "finished_genomes")
    output: os.path.join(working_dir, "finished_referenceseeker")
    benchmark: os.path.join(working_dir, "benchmarks/refseeker.bmk")
    params:
        working_dir = working_dir
    log: os.path.join(working_dir, "logs/run_referenceseeker.log")
    conda:
        'referenceseeker'
    threads: int(config["threads"])
    shell:
        """
        #cd {params.working_dir}
        mkdir -p referenceseeker
        cd referenceseeker
        cp -r ../genomes/genomes/ .
        cd genomes
        for i in `ls *fasta|sed 's/.fasta//g'`;do referenceseeker ~/M5P_databases/referenceseeker/bacteria-refseq $i.fasta > ../$i.RS.tsv --threads {threads};done
        cd ..
        for a in *.tsv; do b=$(basename $a .RS.tsv);  sed -e '1,1d' $a | sed "s/^/\t$b /" > $b.filename.tsv; done
        echo MAG_ID$'\t'ID$'\t'Mash_Distance$'\t'ANI$'\t'Con_DNA$'\t'Taxonomy_ID$'\t'Assembly_Status$'\t'Organism > refseeker.tsv
        for j in *.filename.tsv; do if [ -s "$j" ]; then cat *filename* >> refseeker.tsv; fi; done
        sed -i 's/^\t//g' refseeker.tsv
        rm *filename*
        mv refseeker.tsv ..
        cd ..
        #echo 'Refseeker_complete' > {log}
        touch {output}
        """