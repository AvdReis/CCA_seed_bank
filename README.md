# CCA_seed_bank
This is the code used to analyze the raw fastq files using bash scripts, and subsequently compile/filter the data in R.

## Step 1:
File: Cutadapt.txt
Cutadapt was used to remove the primers.
Gene regions were pooled for indexing: 16S-ITS; 23S-COI; rbcL-tufA; 18S.
A shortcut to making the path.txt needing in next step can be found here.  

## Step 2:
File: P1_Q2.bash
Qiime 2 was used to vizualise reads (using the online Qiime view tool) for each gene region, respectively.

## Step 3
File: P2_Q2.bash
DADA2 within Qiime 2 is used to denoise and remove sequencing errors from the raw data.
Need to adjust trunc values for each gene as per resulting data from P1_Q2.bash
--p-trunc-len-f x \ #16S 270;ITS 270; 23S 270; COI 270; rbcL 270; tufA 250; 18S 270
--p-trunc-len-r x \ #16S 230; ITS 190; 23S 220; COI 230; rbcL 230; tufA 230; 18S 250

## Step 4
File: Blast.bash
Assigning taxonomy using the GenBank nucleotide DNA database using megaBLAST.

## Step 5
File: Step5_CCA_R_code.R
Code used for data filtering, analyses and figure production in R.

## Step 6
File: Step6_tree.txt
Code used for making phylogenetic tree files.
