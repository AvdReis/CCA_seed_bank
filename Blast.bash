#!/bin/bash
#SBATCH -J Blast_CCA
#SBATCH -A xxx #fill in
#SBATCH --time=24:00:00
#SBATCH --mem=100GB
#SBATCH --cpus-per-task=15
#SBATCH --mail-user=xxx #fill in email address for run notifications
#SBATCH --mail-type=ALL

module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2022-07

export TMPDIR=/path/to/temp/directory/can/be/stored/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

QUERIES="./dna-sequences.fasta"
FORMAT="6 qseqid sskingdoms qstart qend qlen qseq qcovs pident sseqid sgi sacc sstart send staxids sscinames stitle length evalue bitscore qcovhsp"

BLASTOPTS="-evalue 0.001 -max_target_seqs 1 -task megablast -perc_identity 0.9"
BLASTAPP="blastn"
DB="nt"

$BLASTAPP $BLASTOPTS -db $DB -query $QUERIES -outfmt "$FORMAT" \
    -out $QUERIES.$DB.$BLASTAPP -num_threads $SLURM_CPUS_PER_TASK
