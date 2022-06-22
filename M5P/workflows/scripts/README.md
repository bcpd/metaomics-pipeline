## Initiation scripts for automation of Atlas-Metagenomes

### concatReads.py

usage: concatReads.py [-h] -i INPUT [-o OUTPUT] [-d DOWNRATE] [-r DOWNREADS]

optional arguments:

  -h, --help            show this help message and exit
  
  -i INPUT, --input INPUT
                        Input directory containing read files
  
  -o OUTPUT, --output OUTPUT
                        Output file basename (optional)
  
  -d DOWNRATE, --downrate DOWNRATE
                        Downsampling rate (optional)
  
  -r DOWNREADS, --downreads DOWNREADS
                        Downsampling read depth (optional)


### formatSamples.py

usage: formatSamples.py [-h] -s SAMPLES [-m METADATA] [-a]

options:
  
  -h, --help            show this help message and exit
  
  -s SAMPLES, --samples SAMPLES
                        Input samples.tsv file path
  
  -m METADATA, --metadata METADATA
                        Metadata file (optional)
  
  -a, --all             Put all samples in same bin group


### run_bakta.sh

usage: run_bakta.sh -i input [-d input] [-o  output directory} [-d databases directory]

Arguments:
    -i                  directory containing MAGs
    -o                  output directory containing per-MAG annotations
    -d                  directory containing annotation databases


### run_bakta.sh

usage: run_bakta.sh -i input [-d input] [-o  output directory} [-d databases directory]

Arguments:
    -i                  directory containing MAGs
    -o                  output directory containing per-MAG annotations
    -d                  directory containing annotation databases


### run_referenceseeker.sh

usage: run_referenceseeker.sh [-i input] [-o  output directory} [-d databases directory]
    -i                  directory containing MAGs
    -o                  output directory containing per-MAG annotations
    -d                  directory containing classification databases
 
 
### create_grist_config_file.py
 
 
usage: pyton create_grist_config_file.py [-s samples_folder] [-d grist database folder] [-o output configuration file ] [-m maximum memory]
    -s                  Folder with metatranscriptomics fastq (path)
    -d                  Directory containing grist databases (path)
    -o                  Config file used to run grist
    -m                  Maximum RAM memory, default=10e9