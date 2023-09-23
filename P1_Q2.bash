#!/bin/bash
#SBATCH -J Q2_P1_CCA
#SBATCH -A xxx #fill in
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=2

module load QIIME2/2022.2

#Importing our files into qiime format
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./path.txt \
  --output-path CCA_pe_demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

#Visualization of sequence quality - sequence quality control - Viewing a summary of joined data with read quality
qiime demux summarize \
  --i-data CCA_pe_demux.qza \
  --o-visualization CCA_pe_demux.qzv \
  --verbose 
