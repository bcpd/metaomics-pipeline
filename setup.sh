# Required programs 


# Get and install miniconda:
wget -L https://repo.anaconda.com/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh
bash  Miniconda3-py37_4.11.0-Linux-x86_64.sh

# Get Snakemake using recommended distributions

conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake

# Get docker
curl -fsSL https://get.docker.com -o get-docker.sh
sh get-docker.sh



# Install conda environments
## These environments are the core of the pipeline. We can run the rules and point to them directly.

# Databases location
## We use a folder in the main home directory to store the databases files used by the different programs.
## This folder is called  M5P_databases 

mkdir ~/M5P_databases

## Atlas
conda install -c bioconda -c conda-forge metagenome-atlas


## Bakta
conda env create -n bakta
conda activate bakta
conda install -c conda-forge -c bioconda bakta
mkdir ~/M5P_databases/bakta
bakta_db download --output ~/M5P_databases/bakta
conda deactivate


## Praxis
git clone https://github.com/davidlevybooth/Praxis.git
cd Praxis
conda env create --name Praxis --file envs/environment.yaml
conda activate Praxis
python setup.py install
conda deactivate
cd ~


## Grist
conda create -y -n grist python=3.9 pip
conda activate grist
python -m pip install genome-grist
mkdir -p ~/M5P_databases/
cd ~/M5P_databases/grist
# Get Grist database GTDB R07-RS207 all genomes (318k)
curl -JLO https://osf.io/k2u8s/download  # 31K-mer database
conda deactivate
cd ~


## Referenceseeker
conda env create -n referenceseeker
conda activate referenceseeker
conda install -c bioconda referenceseeker
mkdir ~/M5P_databases/referenceseeker
wget -L https://zenodo.org/record/4415843/files/bacteria-refseq.tar.gz
tar -xzf *
conda deactivate
cd ~


## M5P enviroment
conda env create --name M5P --file M5P/workflows/envs/environment.yaml