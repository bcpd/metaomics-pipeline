#! /bin/bash
source ./etc/profile.d/conda.sh
wget https://raw.githubusercontent.com/shafferm/DRAM/master/environment.yaml
conda env  create --name DRAM -f environment.yaml
conda activate DRAM
DRAM-setup.py prepare_databases --output_dir DRAM_data --skip_uniref --verbose --threads 10 &> /logs/download_dram.log
DRAM-setup.py update_description_db
DRAM-setup.py export_config --output_file /logs/DRAM_config.log
