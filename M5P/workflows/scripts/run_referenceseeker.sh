#!/bin/bash -e

# Arguments:
# i: directory containing MAGs
# o: output directory containing per-MAG and concatenated referenceseeker annotations
# d: directory containing classification database

while getopts i:o:d: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
        d) db_dir=${OPTARG};;
    esac
done

mkdir -p ${input_dir}
cd ${input_dir}

for f in *.fasta; do
  j=$(basename $f .fasta)
  referenceseeker $(db_dir} ${f} --threads 8 > ${j}.RS.tsv
done

echo MAG_ID$'\t'ID$'\t'Mash_Distance$'\t'ANI$'\t'Con_DNA$'\t'Taxonomy_ID$'\t'Assembly_Status$'\t'Organism > refseeker.tsv

for a in *.tsv; do
  b=$(basename $a .RS.tsv)
  sed -e '1,1d' ${a} | sed "s/^/\t${b} /" > ${b}.filename.tsv
done

for j in *.filename.tsv; do
  if [ -s "$j" ]; then
        cat *filename* >> refseeker.tsv
  fi
done

rm *filename*

mkdir -p ${output_dir}
mv *RS.tsv ${output_dir}
mv refseeker.tsv ${output_dir}
