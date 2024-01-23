#!/bin/bash

cd /home/jmarin/work/Hemato/pip2
mkdir stats
mkdir whamg


# Create the files list_sam.txt in the directory
ls mapping_2/*.sam > list.sam.txt

conda activate samtools-1.14
while read line
do
	filename=$line
	filename="${filename//.aln-pe.sam}"
	
	samtools view -bS $line > ${filename}.bam
	
	# name order to coordinate order
	samtools sort -O bam -o ${filename}.bam -T tmp ${filename}.bam
	
	# make index
	samtools index -b ${filename}.bam

done < list.sam.txt
conda deactivate	
	
# 
# 	# clean up read pairing information and flags
# 	samtools fixmate -O bam $line ${filename}.aln-fixmate.bam
# 
# 	# name order to coordinate order
# 	samtools sort -O bam -o ${filename}.aln-fixmate.bam -T tmp ${filename}.aln-fixmate.bam
# 	
# 	# get statistics about alignments
# 	outdir="${filename/mapping_2/stats}"
# 	samtools index -b ${filename}.aln-fixmate.bam
# 	samtools stats ${filename}.aln-fixmate.bam > ${outdir}.stats_summary.txt
# 	samtools flagstat ${filename}.aln-fixmate.bam > ${outdir}.stats.txt
# 
# 	sequence="${outdir//stats\/}"
# 	sequence="${sequence//_L001_R1_001_val_1.fq.gz}"
# 	
# 	samtools depth ${filename}.aln-fixmate.bam > ${outdir}.depth.txt
# 	
# 	cov=`samtools depth -a ${filename}.aln-fixmate.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'`
# 	echo $sequence $cov >> stats/CoverageBreadth.txt
# 	
#     dep=`samtools depth -a ${filename}.aln-fixmate.bam | awk '{c++;s+=$3}END{print s/c}'`
#     echo $sequence $dep >> stats/MeanReadDepth.txt
#     
# 	nb=`grep ^"SN	raw total sequences:" ${outdir}.stats_summary.txt | cut -d ":" -f2`
# 	echo $sequence $nb >> stats/NbSeq.txt
# 	
# 	qual=`grep ^"SN	average quality:" ${outdir}.stats_summary.txt | cut -d ":" -f2`
# 	echo $sequence $qual >> stats/Quality.txt
# 
# done < list.sam.txt
# conda deactivate


#### General structural variant (SV) discovery

# Create the files list_bam.txt in the directory
ls mapping_2/*.bam > list.bam.txt

conda activate wham-1.8.0.1

export EXCLUDE="caca"
whamg -x 2 -e $EXCLUDE -a ../snippy/ref/patient2.fasta -f list.bam.txt > whamg/patient2.vcf  2> whamg/patient2.err

conda deactivate













