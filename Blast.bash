#!/bin/bash
#SBATCH -J LU_2022
#SBATCH -A uoa03362
#SBATCH --time=24:00:00
#SBATCH --mem=100GB
#SBATCH --cpus-per-task=15
#SBATCH --mail-user=avan398@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2022-07

export TMPDIR=/nesi/nobackup/uoa03362/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

QUERIES="./dna-sequences.fasta"
FORMAT="6 qseqid sskingdoms qstart qend qlen qseq qcovs pident sseqid sgi sacc sstart send staxids sscinames stitle length evalue bitscore qcovhsp"

BLASTOPTS="-evalue 0.001 -max_target_seqs 1 -task megablast -perc_identity 0.9"
BLASTAPP="blastn"
DB="nt"

$BLASTAPP $BLASTOPTS -db $DB -query $QUERIES -outfmt "$FORMAT" \
    -out $QUERIES.$DB.$BLASTAPP -num_threads $SLURM_CPUS_PER_TASK
