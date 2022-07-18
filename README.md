# M5P: **MBI Microbial Mixtures Metagenomics and Metatranscriptomics pipeline**

Snakemake-based workflow for the taxonomic and functional annotation of metagenomics and metatranscriptomics datasets.

## Dependencies:
- Docker
- Miniconda
- [atlas-metagenome](https://github.com/metagenome-atlas/atlas)
- [Praxis](https://github.com/davidlevybooth/Praxis)
- [DRAM](https://github.com/shafferm/DRAM)
- [Bakta](https://github.com/oschwengers/bakta)
- [Referenceseeker](https://github.com/oschwengers/referenceseeker)
- [sourmash](https://github.com/sourmash-bio/sourmash)
- [genome-grist](https://github.com/dib-lab/genome-grist)

## Installation:
- Clone the M5P repository, ie: `gh repo clone MicrobiomeInsights/metaomics-pipeline`
- Alternatively you can download the code and unzip it in an installation folder using `unzip metaomics-pipeline'
- Change into the repo: `cd metaomics-pipeline`
- Install the M5P CLI: `python setup.py install`
- Run `bash setup.sh` to create environments and download databases

## Resource requirements:
- OS: Linux (Ubuntu >= 14.04 LTS).
- Disk space: M5P requires at least 200GB of disk space for databases and sequencing files.
- Memory: Both metagenomics and metatranscriptomics analyses can be done with 128 GB of RAM or less.
- Time:
  - *Installation and setup*: several hours.
  - *Analysis*: largely dependent on the size of the metagenome, number of samples, and anlysis mode, but expect a few hours for the metagenomics-only mode. The processing of multiple data sets can be done in parallel (see below).

## Instructions

The first step is to activate the M5P conda environments to have access to the required libraries (`conda activate M5P`). Then, decide if you want to run the metagenomics pipeline or the metatranscriptomics one. 

You can run the pipeline by specifying the different parameters or by using an editable configuration file e.g. `M5P -c M5P_config.yaml`. The configuration file is included in the repository. A template that you can edit can be created using the command `M5P -a`.
We suggest that you name the configuration file as *M5P_config.yaml*. Do not name it as *config.yaml* because a file with that name is created in the process for one of the components.


**Example use for metagenomics**

`M5P -e metagenomics - w . -i input_folder -t 8 -j 6 -M 128000000 -m exampleMetadata.tsv -e metagenomics`

or 

`M5P -c M5P_config.yaml`


**Example use for metatranscriptomics**

`M5P -e metatranscriptomics -w . -x input_folder -t 8 -j 8 -M 128000000`

Note: for running differential abundance tests with DESeq2, please enable the relevant options using the M5P_config.yaml file.



## Usage:
```
M5P [-h] [-w WORKING_DIR] [-i FASTQ_METAGENOMICS] [-x FASTQ_METATRANSCRIPTOMICS] [-r MERGED_READS] [-m METADATA_PATH] [-t THREADS]
           [-p MERGED_PREFIX] [-b BIN_ALL] [-j JOBS] [-c CONFIGFILE] [-a GETCONFIG] 

Pipeline configuration options:

optional arguments:
  -h, --help            show this help message and exit
  -w WORKING_DIR, --working_dir WORKING_DIR
                        working directory (path)
  -i FASTQ_METAGENOMICS, --fastq_metagenomics FASTQ_METAGENOMICS
                        Directory containing metagenomics fastq files (path)
  -x FASTQ_METATRANSCRIPTOMICS, --fastq_metatranscriptomics FASTQ_METATRANSCRIPTOMICS
                        Directory containing metatranscriptomics fastq files
                        (path)
  -r MERGED_READS, --merged_reads MERGED_READS
                        Merge reads for co-assembly (True or False)
  -m METADATA_PATH, --metadata_path METADATA_PATH
                        Metadata file (path)
  -t THREADS, --threads THREADS
                        The number of threads the pipeline is allowed to use
                        (integer)
  -M MAX_MEMORY, --max_memory MAX_MEMORY
                        The amount of memory provided to the assembler. Enter
                        in byte format.
  -e EXPERIMENT_TYPE, --experiment_type EXPERIMENT_TYPE
                        Either metagenomics (default), metatranscriptomics, or
                        both
  -k EXPERIMENTAL_CONTRAST, --experimental_contrast EXPERIMENTAL_CONTRAST
                        Experimental contrast as used in R formula
  -g EXPERIMENTAL_DESIGN, --experimental_design EXPERIMENTAL_DESIGN
                        Experimental design as used in R formula
  -p MERGED_PREFIX, --merged_prefix MERGED_PREFIX
                        Prefix for merged reads files (string)
  -b BIN_ALL, --bin_all BIN_ALL
                        Put all reads files in same bin group (bool)
  -j JOBS, --jobs JOBS  number of jobs
  -c CONFIGFILE, --configfile CONFIGFILE
                        optional yaml file containing all of the above
                        configuration details.
  -a, --GETCONFIG       Prints a yaml template with the default configuration
   ```
