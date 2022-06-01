rule get_dram:
    shell:
        "docker pull continuumio/miniconda3"
        " "
        "docker run --name DRAM continuumio/miniconda3"
        " "
        "docker exec DRAM mkdir -p data out scripts logs"
        " "
        "docker cp ../scripts/DRAM_setup.sh DRAM:/scripts"
        " "
        "docker exec DRAM /bin/bash /scripts/DRAM_setup.sh"

rule annotate_genomes:
    shell:
        "docker cp /MAGs/fasta/*fasta DRAM:/genomes/"
        " "
        "docker cp ../scripts/annotate_genomes.sh DRAM:/scripts"
        " "
        "docker exec -t DRAM /bin/bash /scripts/annotate_genomes.sh"


rule copy_DRAM_annotations:
   shell:
       "docker cp DRAM:out/ MAGS/functional_annotations/"
       " "
       "docker cp DRAM:logs/* /logs/"
