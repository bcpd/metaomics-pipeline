# Required programs


# Get and install miniconda:
wget -L https://repo.anaconda.com/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh
bash  Miniconda3-py37_4.11.0-Linux-x86_64.sh

# Get Snakemake using recommended distributions

conda install -y -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake

# Get docker
curl -fsSL https://get.docker.com -o get-docker.sh
sh get-docker.sh


# Databases location
## We use a folder in the main home directory to store the databases files used by the different programs.
## This folder is called  M5P_databases

mkdir ~/M5P_databases


# Install conda environments
## These environments are the core of the pipeline. We can run the rules and point to them directly.

## Atlas
conda create -n atlas
conda activate atlas
mamba install -y -c bioconda -c conda-forge metagenome-atlas=2.9
conda deactivate

## GTDBTK
## Atlas uses gtdbtk for classification of MAGs but sometimes it runs into a problem downloading it
## So we download it directly and create a file that is later used to see if the db is present
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz -P ~/M5P_databases/GTDB_V06/ --no-check-certificate
tar xvzf ~/M5P_databases/GTDB_V06/gtdbtk_r202_data.tar.gz -C ~/M5P_databases/GTDB_V06/ --strip 1
rm ~/M5P_databases/GTDB_V06/gtdbtk_r202_data.tar.gz
touch ~/M5P_databases/GTDB_V06/downloaded_success

## Bakta
conda create -n bakta
conda activate bakta
conda install -y -c conda-forge -c bioconda bakta=1.4.2
mkdir ~/M5P_databases/bakta
bakta_db download --output ~/M5P_databases/bakta
conda deactivate

## Praxis
git clone https://github.com/davidlevybooth/Praxis.git
cd Praxis
mamba env create --name Praxis --file envs/environment.yaml
conda activate Praxis
python setup.py install
conda deactivate
cd ~

## Sourmash
pip install sourmash==4.3.0

## Grist
conda create -y -n grist python=3.9 pip
conda activate grist
python -m pip install genome-grist
mkdir ~/M5P_databases/grist
cd ~/M5P_databases/grist
# Get Grist database GTDB R07-RS207 all genomes (318k)
curl -JLO https://osf.io/k2u8s/download  # 31K-mer database
conda deactivate
cd ~

## Referenceseeker
conda create -n referenceseeker
conda activate referenceseeker
conda install -y -c bioconda referenceseeker=1.8.0
mkdir ~/M5P_databases/referenceseeker
cd ~/M5P_databases/referenceseeker
wget -L https://zenodo.org/record/4415843/files/bacteria-refseq.tar.gz
tar -xzf *
conda deactivate
cd ~

## DRAM
## To correctly setup DRAM you need to have access to three scripts
## that are located in the "M5P/workflows/scripts" folder
## After the installation and setup we create an small empty file in the
## M5P_databases folder to signal to the pipeline that he installation worked

docker pull continuumio/miniconda3
docker run -i -d --name DRAM continuumio/miniconda3
docker exec DRAM mkdir -p data out scripts logs genomes proteins
docker cp M5P/workflows/scripts/DRAM_setup.sh DRAM:/scripts
docker cp M5P/workflows/scripts/DRAM_annotate_genomes.sh DRAM:/scripts
docker cp M5P/workflows/scripts/DRAM_annotate_proteins.sh DRAM:/scripts
docker exec DRAM /bin/bash /scripts/DRAM_setup.sh
touch ~/M5P_databases/DRAM_installed

## dRep enviroment
mamba env create --name drep --file M5P/workflows/envs/drep.yaml
conda activate drep
# After installation, we recommend testing the access to the different required programs using drep check_dependecies

mamba install -y -c bioconda checkm-genome
wget -L https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar
tar -xvf mash-Linux64-v2.3.tar
cp mash-Linux64-v2.3/mash $CONDA_PREFIX/bin
wget -L https://ani.jgi.doe.gov/download_files/ANIcalculator_v1.tgz
tar -xvf ANIcalculator_v1.tgz
cp ANIcalculator_v1/ANIcalculator $CONDA_PREFIX/bin
cp ANIcalculator_v1/nsimscan $CONDA_PREFIX/bin
conda deactivate

## M5P enviroment
mamba env create --name M5P --file M5P/workflows/envs/environment.yaml
# conda activate M5P
# Get zip file with m5p code
# unzip zip file
# Enter folder `cd metaomics-pipeline-main`
# Install: `python setup.py install`
# conda deactivate


