#! /bin/bash

### usage ###
# sh chrom_fasta_rmnewlines_Nov072023_hj.sh [gene] [target chrom] [header line for target chrom + 1] [header line for chrom right after target chrom - 1] 

gene=`echo $1`
chrom=`echo $2`
start=`echo $3`
end=`echo $4`

mkdir -p ${gene}

sed -n "${start},${end}p" /Volumes/hjdrive/thy_n25/thy_recurrent/20231107/spliceai/human_g1k_v37_decoy.fasta > ./${gene}/human_g1k_v37_chr${chrom}.fasta

# remove line
cat ./${gene}/human_g1k_v37_chr${chrom}.fasta | tr -d '\n' > ./${gene}/human_g1k_v37_chr${chrom}_rmnewline.fasta

rm ./${gene}/human_g1k_v37_chr${chrom}.fasta # remove tmp file 
