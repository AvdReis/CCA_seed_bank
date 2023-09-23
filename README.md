# CCA_prelim
This is the code used to analyze the raw fastq files using bash scripts, and subsequently compile/filter the data in R.

## Step 1:
Cutadapt was used to remove the primers.
Gene regions were pooled for indexing:
16S-ITS
23S-COI
rbcL-tufA
18S

## Step 2:
Qiime 2 was used to vizualise reads (using the online Qiime view tool) for each gene region, respectively.
--p-trunc-len-f x \ #16S 270;ITS 270; 23S 270; COI 270; rbcL 270; tufA 250; 18S 270
--p-trunc-len-r x \ #16S 230; ITS 190; 23S 220; COI 230; rbcL 230; tufA 230; 18S 250
