#! /bin/bash
source ./etc/profile.d/conda.sh
conda env  create --name DRAM -f /scripts/dram.yaml
conda activate DRAM
cp /scripts/DRAM-setup.py /opt/conda/envs/DRAM/bin/DRAM-setup.py  # Copy corrected DRAM setup file
DRAM-setup.py prepare_databases --output_dir DRAM_data --skip_uniref --verbose --threads 10 &> /logs/download_dram.log
DRAM-setup.py update_description_db
DRAM-setup.py export_config --output_file /logs/DRAM_config.log
