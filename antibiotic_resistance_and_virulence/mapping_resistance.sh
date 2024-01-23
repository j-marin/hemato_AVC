#!/bin/bash

cd /home/jmarin/work/Hemato/pip2
mkdir mapping_resistance

# Create the files list_sample.txt, list_R1.txt and list_R2.txt in the directory
find trim_results_all94 -maxdepth 1 -name "*R1_001_val_1.fq*" > list_val_R1.txt
sed 's/_R1_001_val_1.fq/_R2_001_val_2.fq/g' list_val_R1.txt > list_val_R2.txt

## (3) Map SNP
conda activate bwa-0.7.17
count=1
while read lineA
    do 
        filename=$lineA
        filename="${filename//_L001_R1_001_val_1.fq.gz}"
        filename="${filename/trim_results_all94/mapping_resistance}"
        
        lineB=`sed -n "$count"p list_val_R2.txt`
        count=`expr $count + 1`
        
 		id=`basename $filename`
 
 		echo $id >> list_id.txt
		
        bwa index ../snippy/ref/RefSeq_resistance.fasta
        bwa mem ../snippy/ref/RefSeq_resistance.fasta -t 3 -R "@RG\tID:$id\tSM:$id" $lineA $lineB > ${filename}.aln-pe.sam
       
done < list_val_R1.txt
conda deactivate

echo 'step 2 complete' > step2.txt


