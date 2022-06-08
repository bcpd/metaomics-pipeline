#!/bin/bash

# adds MAG identifier to each bakta annotation file (MAG{number}) and concatenates annotation files into single annotation file
# Arguments: folder containing bakta annotations. The annotations file is saved in the same folder as the MAG* files.
#!/bin/bash -e

# Arguments:
# i: directory containing MAGs
# o: output directory containing per-MAG annotations
# d: directory containing annotation databases

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

for i in *.fasta; do
  bakta --db ${db_dir} ${i} --output ${output_dir} --threads 8
done

cd ${output_dir}

echo SRA$'\t'Sequence_Id$'\t'Type$'\t'start_position$'\t'end_position$'\t'Strand$'\t'Locus_Tag$'\t'Gene$'\t'Product$'\t'Product$'\t'DbXrefs > bakta.tsv

for f in *.tsv; do
        if [[ $f != *"hypotheticals"* ]] ; then
          i=$(basename $f .tsv)
          sed -e '1,3d' ${f} | sed "s/^/\t${i} /" > ${i}.filename.tsv
          cat *filename* >> bakta.tsv
      fi
done

#mv bakta.tsv ${output_dir}

rm bakta.filename.tsv
