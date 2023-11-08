![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/486d2656-b137-4025-82b9-274564aa42d8)![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/4009b1b1-576d-4670-8139-32113b35143a)# spliceAI_rawscore  
this repository is for showing my method to use raw score of spliceAI.  

# Background  
When we have to predict the consequence of splicing site mutation, we can utilize spliceAI.  

# Biological question  
![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/8796e047-6b15-4286-8954-34fb6417e476)

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



![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/cf560bfe-4497-4883-b69f-3b95cadb3b8a)




![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/33a3da56-3857-4196-96b5-06342e1b628f)


![image](https://github.com/hojinlee-98/spliceAI_rawscore/assets/121307215/7948f152-cf68-4ef9-a0c6-8f408bf57a60)
