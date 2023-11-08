#! /bin/bash

gene=`echo $1`

echo "start" > ./${gene}/gencode_v24lift37_${gene}_hj.txt
grep "\s${gene}\s" 20230708_gencode_v24lift37_hj.txt | awk -F '\t' '{print $5}' | sort -n >> ./${gene}/gencode_v24lift37_${gene}_hj.txt
echo "end" >> ./${gene}/gencode_v24lift37_${gene}_hj.txt
grep "\s${gene}\s" 20230708_gencode_v24lift37_hj.txt | awk -F '\t' '{print $6}' | sort -n -r >> ./${gene}/gencode_v24lift37_${gene}_hj.txt
