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
