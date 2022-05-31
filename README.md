# metaomics-pipeline

**MBI Microbial Mixtures Metagenomics and Metatranscriptomics pipeline**

Snakemake-based workflow for the taxonomic and functional annotation of metagenomics and metatranscriptomics datasets.

Dependencies:
- [atlas-metagenome](https://github.com/metagenome-atlas/atlas)
- [Praxis](https://github.com/davidlevybooth/Praxis)
- [DRAM](https://github.com/shafferm/DRAM)

Installation:
- Run `bash setup.sh` to create environments and download databases
- Clone the M5P repo, ie: `gh repo clone MicrobiomeInsights/metaomics-pipeline`
- Change into the repo: `cd metaomics-pipeline`
- Install the M5P CLI: `python setup.py install`

usage: M5P [-h] [-w WORKING_DIR] [-i FASTQ_DIR] [-d DATABASE_DIR] [-r MERGED_READS] [-m METADATA_PATH] [-t THREADS]
           [-p MERGED_PREFIX] [-b BIN_ALL] [-c CONFIGFILE]

pipeline configuration

options:
  -h, --help            show this help message and exit

  -w WORKING_DIR, --working_dir WORKING_DIR
                        working directory (path)

  -i FASTQ_DIR, --fastq_dir FASTQ_DIR
                        Directory containing fastq files (path)

  -d DATABASE_DIR, --database_dir DATABASE_DIR
                        Directory containing ATLAS databases (path)

  -r MERGED_READS, --merged_reads MERGED_READS
                        Merge reads for co-assembly (True or False)

  -m METADATA_PATH, --metadata_path METADATA_PATH
                        Metadata file (path)

  -t THREADS, --threads THREADS
                        The number of threads the pipeline is allowed to use (integer)

  -p MERGED_PREFIX, --merged_prefix MERGED_PREFIX
                        Prefix for merged reads files (string)

  -b BIN_ALL, --bin_all BIN_ALL
                        Put all reads files in same bin group (bool)

  -c CONFIGFILE, --configfile CONFIGFILE
                        optional yaml file containing all of the above configuration details.

