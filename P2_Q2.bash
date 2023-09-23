#!/bin/bash
#SBATCH -J Q2_P2_16S_CCA
#SBATCH -A uoa03362
#SBATCH --time=10:00:00
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=8

module load QIIME2/2022.2

#this may not be necessary in all cases
export TMPDIR=/path/to/temp/directory/can/be/stored/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

#Quality and chimera filtering using DADA2
## No need to trim as already removed primers using cutadapt
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./CCA_pe_demux.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 270 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 230 \
  --p-chimera-method consensus \
  --p-min-fold-parent-over-abundance 1 \
  --p-n-threads 0 \
  --o-representative-sequences rep_seqs_dada2.qza \
  --o-table table_dada2.qza \
  --o-denoising-stats CCA-stats-dada2.qza \
  --verbose

#Looking at the feature table and feature data summary
qiime feature-table summarize \
  --i-table table_dada2.qza \
  --o-visualization table_dada2.qzv

qiime feature-table tabulate-seqs \
  --i-data rep_seqs_dada2.qza \
  --o-visualization rep_seqs_dada2.qzv

#Had issues with this once and have left it in since
#if no features/OTUs are present in the table this creates a problem trying to read in the table in R
qiime feature-table filter-samples \
 --i-table table_dada2.qza \
 --p-min-features 1 \
 --o-filtered-table rm_0_features_table.qza

# Exporting table
qiime tools export \
  --input-path ./rm_0_features_table.qza\
  --output-path ./

#Exporting rep seqs
qiime tools export --input-path ./rep_seqs_dada2.qza --output-path ./

#looking at denoising stats from dada2
  qiime metadata tabulate \
  --m-input-file CCA-stats-dada2.qza \
  --o-visualization CCA-stats-dada2.qzv
