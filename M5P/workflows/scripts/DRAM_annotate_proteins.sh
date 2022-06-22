#! /bin/bash
source ./etc/profile.d/conda.sh
conda activate DRAM

# Run dram annotate
DRAM.py annotate_genes -i '/genomes/*faa'  -o /out/annotations/ --threads 10 --verbose &> /logs/dram_protein_annotation.log

# Run dram distill
DRAM.py distill --input_file /out/annotations/annotations.tsv --output_dir /out/annotations/distill  &> /logs/dram_protein_distil.log
