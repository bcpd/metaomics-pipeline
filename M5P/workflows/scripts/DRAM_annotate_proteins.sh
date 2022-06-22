#! /bin/bash
source ./etc/profile.d/conda.sh
conda activate DRAM

# Run dram annotate
DRAM.py annotate_genes -i '/genomes/*faa'  -o /out/annotations/ --threads 10 --verbose &> /logs/dram_protein_annotation.log

# Run dram distill
if [ -f /out/annotations/rrnas.tsv ] && [ -f /out/annotations/trnas.tsv ]
then
   DRAM.py distill --input_file /out/annotations/annotations.tsv --rrna_path /out/annotations/rrnas.tsv --trna_path /out/annotations/trnas.tsv --output_dir /out/annotations/distill  &> /logs/dram_distil.log
elif [ -f /out/annotations/rrnas.tsv ]
then
   DRAM.py distill --input_file /out/annotations/annotations.tsv --rrna_path /out/annotations/rrnas.tsv --output_dir /out/annotations/distill  &> /logs/dram_distil.log
else
   DRAM.py distill --input_file /out/annotations/annotations.tsv --trna_path /out/annotations/trnas.tsv --output_dir /out/annotations/distill  &> /logs/dram_distil.log
fi
