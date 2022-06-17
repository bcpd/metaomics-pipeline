# Get miniconda:
# Note sure which one
wget -L https://repo.anaconda.com/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh
#wget -L https://repo.anaconda.com/miniconda/Miniconda3-py38_4.11.0-Linux-x86_64.sh
#wget -L https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh

bash  Miniconda3-py37_4.11.0-Linux-x86_64.sh
#bash Miniconda3-py38_4.11.0-Linux-x86_64.sh
#bash Miniconda3-py39_4.11.0-Linux-x86_64.sh

# Get Snakemake using recommended distributions

conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake

# Get docker

curl -fsSL https://get.docker.com -o get-docker.sh
sh get-docker.sh


# create folder structures



# Install conda environments
# THse environments are the core of the pipeline. We can run the rules and point to them directly.
# The alternative is to create enviroment files (yml) so conda can create them as needed, though these environments are not automatically purged after use.

## Atlas
conda install -c bioconda -c conda-forge metagenome-atlas

## Bakta
conda install -c conda-forge -c bioconda bakta
## download bakta db
bakta_db download --output {database_dir}

## referenceseeker
conda install -c bioconda referenceseeker

mkdir -p {database_dir}
cd {database_dir}
wget https://zenodo.org/record/4415843/files/bacteria-refseq.tar.gz
tar -xzf bacteria-refseq.tar.gz

## Praxis
git clone https://github.com/davidlevybooth/Praxis.git
cd Praxis
conda env create --name Praxis --file envs/environment.yaml
conda activate Praxis
python setup.py install
cd ..


## Grist
conda create -y -n grist python=3.9 pip
conda activate grist
python -m pip install genome-grist
curl -JLO https://osf.io/9gpck/download
curl -JLO https://osf.io/k2u8s/download
curl -JLO https://osf.io/ubt7p/download

## M5P enviroment
conda env create --name M5P --file M5P/workflows/envs/environment.yaml
