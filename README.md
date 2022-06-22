# metaomics-pipeline

**MBI Microbial Mixtures Metagenomics and Metatranscriptomics pipeline**

Snakemake-based workflow for the taxonomic and functional annotation of metagenomics and metatranscriptomics datasets.

Dependencies:
- Docker
- Miniconda
- [atlas-metagenome](https://github.com/metagenome-atlas/atlas)
- [Praxis](https://github.com/davidlevybooth/Praxis)
- [DRAM](https://github.com/shafferm/DRAM)
- [Bakta](https://github.com/oschwengers/bakta)
- [Referenceseeker](https://github.com/oschwengers/referenceseeker)


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

  -i FASTQ_METAGENOMICS, --fastq_metagenomics FASTQ
                        Directory containing metagenomics fastq files (path)

  -x FASTQ_METATRANSCRIPTOMICS, --fastq_metatranscriptomics FASTQ
                        Directory containing metatranscriptomics fastq files (path)

  -r MERGED_READS, --merged_reads MERGED_READS
                        Merge reads for co-assembly (True or False), default=True

  -m METADATA_PATH, --metadata_path METADATA_PATH
                        Metadata file (path)

  -t THREADS, --threads THREADS
                        The number of threads the pipeline is allowed to use (integer), default=2

  -M MAXIMUM_MEMORY, --max_memory MAXIMUM_MEMORY
                        The amount of memory provided to the assembler. Enter in byte format


  -p MERGED_PREFIX, --merged_prefix MERGED_PREFIX
                        Prefix for merged reads files (string)

  -b BIN_ALL, --bin_all BIN_ALL
                        Put all reads files in same bin group (bool)

  -c CONFIGFILE, --configfile CONFIGFILE
                        optional yaml file containing all of the above configuration details. Default="M5P_config.yaml"

  -j JOBS, --jobs JOBS
                        Number of jobs to run, default=2
