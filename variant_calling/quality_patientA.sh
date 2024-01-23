#!/bin/bash

cd /home/jmarin/work/Hemato/pip2

# Create the files list_sample.txt, list_R1.txt and list_R2.txt in the directory
ls /home/jmarin/save/Hemato/fastq_clean/patient2/*.fastq* > list.samples.txt
find /home/jmarin/save/Hemato/fastq_clean/patient2 -maxdepth 1 -name "*R1_001.fastq*" > list_R1.txt
sed 's/_R1_001.fastq/_R2_001.fastq/g' list_R1.txt > list_R2.txt

## (1) Quality
# Run fastqc
mkdir fastqc_results_2

conda activate fastqc-0.11.9
while read line
do
	fastqc -o fastqc_results_2 $line
done < list.samples.txt
conda deactivate

# Run multiqc
mkdir multiqc_results_2
cd multiqc_results_2

conda activate multiqc-1.11
multiqc -ip ../fastqc_results_2/
conda deactivate

## (2) Trimming by pairs of reads
cd ..
mkdir trim_results_2

conda activate trim-galore-0.6.8
count=1
while read lineA
    do 
        lineB=`sed -n "$count"p list_R2.txt`
        count=`expr $count + 1`
        trim_galore -q 30 --illumina --paired -o trim_results_2 --length 50 $lineA $lineB

done < list_R1.txt
conda deactivate

echo 'step 1 complete' > step1.txt

