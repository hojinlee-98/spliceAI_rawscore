# spliceAI_rawscore  
this repository is for showing my method to use raw score of spliceAI.  

# Background  
When we have to predict the consequence of splicing site mutation, we can utilize spliceAI.  

# Biological question  
![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/b8fc9d66-ce84-4652-95d6-74cbc0857264)


# Method  
```
grep -n '>' /Volumes/hjdrive/thy_n25/thy_recurrent/20231107/spliceai/human_g1k_v37_decoy.fasta > 20231107_chrom_header_line.txt

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

```
![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/ec5e5670-d52c-49e4-a106-636ff9f330aa)


```
gene=`echo $1`
echo "start" > ./${gene}/gencode_v24lift37_${gene}_hj.txt
grep "${gene}" 20230708_gencode_v24lift37_hj.txt | awk -F '\t' '{print $5}' | sort -n >> ./${gene}/gencode_v24lift37_${gene}_hj.txt
echo "end" >> ./${gene}/gencode_v24lift37_${gene}_hj.txt
grep "${gene}" 20230708_gencode_v24lift37_hj.txt | awk -F '\t' '{print $6}' | sort -n -r >> ./${gene}/gencode_v24lift37_${gene}_hj.tx![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/5da9b5c7-dc82-4045-b815-97e010766645)

(spliceai) hojin@ihojin-ui-iMac PRIM2 % cat gencode_v24lift37_PRIM2_hj.txt 
start
57179602
57182421
57182421
end
57513375
57191168
57191161

python spliceai_prediction_wt_mt_Nov072023_hj.py 6 PRIM2 57372356 G C 57179602 57513375![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/a72f8a66-4592-4a96-aab6-23a2e1ea9e64)

```
![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/6ba2d37d-d5ec-4184-bcd6-fe646bfbd5b4)

![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/e932012b-9cb9-4972-bf15-22599ead8f5e)


# Conclusion
The variants located in splicing site remove exon 8. 
