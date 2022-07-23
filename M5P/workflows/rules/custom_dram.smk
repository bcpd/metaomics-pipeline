#!/usr/bin/python3


rule annotate_genomes:
    '''
    Run DRAM inside the conda image, uses the predicted proteins as input
    '''
    input:
       atlas_genome_complete = os.path.join(working_dir, "finished_genomes")
    output: os.path.join(working_dir, "finished_DRAM_annotate")
    params:
        script = os.path.join(config["parent_dir"], "workflows/scripts/DRAM_annotate_proteins.sh")
    log: os.path.join(working_dir, "logs/DRAM_annotate.log")
    shell:
        """
        docker start DRAM || true
        docker exec DRAM rm -fr /genomes/ ||true
        docker exec DRAM mkdir /genomes/ ||true
        docker exec DRAM rm -fr /out/ ||true
        docker exec DRAM mkdir -p /out/ ||true
        docker exec DRAM rm -fr /logs/ ||true
        docker exec DRAM mkdir -p /logs/ ||true
        for i in genomes/annotations/genes/MAG*faa; do docker cp $i DRAM:/genomes;done
        docker cp {params.script} DRAM:/scripts  || true
        docker exec -t DRAM /bin/bash /scripts/DRAM_annotate_proteins.sh || true
        touch {output}
        touch {log}
        docker stop DRAM ||true
        """


rule copy_DRAM_annotations:
    '''
    Copies the DRAM output files from the docker image to the local directory
    '''
    input: os.path.join(working_dir, "finished_DRAM_annotate")
    output: os.path.join(working_dir, "finished_DRAM")
    params:
        output_folder= os.path.join(working_dir, "DRAM")
    log: os.path.join(working_dir, "logs/DRAM_copy_results.log")
    shell:
        """
        docker start DRAM || true
        mkdir -p {params.output_folder}
        docker cp DRAM:out/annotations {params.output_folder} ||true
        docker cp DRAM:logs/* /logs/ ||true
        docker stop DRAM     || true
        touch {output}
        """


rule annotate_genomes2:
    '''
    Run DRAM inside the conda image, uses genomes generared by grist as input
    '''
    input:
        finished_grist = os.path.join(working_dir, "finished_grist"),
        grist_cleaned =  os.path.join(working_dir, "cleaned_after_grist")
    output: os.path.join(working_dir, "finished_DRAM_annotate_reference_genomes")
    params:
        script = os.path.join(config["parent_dir"], "workflows/scripts/DRAM_annotate_genomes.sh"),
        working_dir = working_dir,
    log: "logs/DRAM_annotate.log"
    shell:
        """
        cd {params.working_dir}
        docker start DRAM || true
        docker exec DRAM rm -fr /genomes/ ||true
        docker exec DRAM mkdir /genomes/ ||true
        docker exec DRAM rm -fr /out/ ||true
        docker exec DRAM mkdir -p /out/ ||true
        docker exec DRAM rm -fr /logs/ ||true
        docker exec DRAM mkdir -p /logs/ ||true
        for i in reference_genomes/*fna; do docker cp $i DRAM:/genomes;done
        docker cp {params.script} DRAM:/scripts  || true
        docker exec -t DRAM /bin/bash /scripts/DRAM_annotate_genomes.sh || true
        docker cp DRAM:out/annotations reference_genomes ||true
        docker cp DRAM:logs/* /logs/ ||true
        touch {log}
        cd -
        touch {output}
        docker stop DRAM ||true
        """