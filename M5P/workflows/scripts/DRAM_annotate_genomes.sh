#! /bin/bash
source ./etc/profile.d/conda.sh
conda activate DRAM

# Run dram annotate
DRAM.py annotate -i '/genomes/*fasta'  -o /out/annotations/ --threads 10 --min_contig_size 1000 --verbose &> /logs/dram_annotation.log

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

# This part did not work
#cd /out/annotations/
#python /scripts/DRAM_get_all_modules_mod.py annotations.tsv kegg_modules.tsv /logs/get_all_modules.log
