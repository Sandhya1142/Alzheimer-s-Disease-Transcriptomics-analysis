#!/bin/bash

#SBATCH --job-name=fastq_download
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --output=fastq_url_download_%j.log

module load anaconda/2024.10
source /iitgn/apps/packages/anaconda/2024.10/etc/profile.d/conda.sh
conda activate rnaseq_env


set -e 

mkdir -p /iitgn/homedirs/sandhya/Rush_AD/URLS_fastq

echo "URL Download in process🌻"
wget -i /iitgn/homedirs/sandhya/Rush_AD/data_processing/fastq.URLs \
     -P	/iitgn/homedirs/sandhya/Rush_AD/URLS_fastq \
     -- continue

echo "Download done🌻"



