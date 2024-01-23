#!/bin/bash

cd /home/jmarin/work/Hemato/pip2
mkdir Fasta
mkdir stats

# Create the files list_sam.txt in the directory
ls mapping_resistance/*.sam > list.sam.txt

conda activate samtools-1.14
while read line
do
	filename=$line
	filename="${filename//.aln-pe.sam}"

	# keep only uniquely mapped reads
	samtools view -b -F 0x0100 -F 0x0400 -F 0x800 $line > ${filename}.unique.sam
	
	# clean up read pairing information and flags
	samtools fixmate -O bam ${filename}.unique.sam ${filename}.aln-fixmate.bam

	# name order to coordinate order
	samtools sort -O bam -o ${filename}.aln-fixmate.bam -T tmp ${filename}.aln-fixmate.bam
	
	# get statistics about alignments
	outdir="${filename/mapping_resistance/stats}"
	samtools index -b ${filename}.aln-fixmate.bam
	samtools stats ${filename}.aln-fixmate.bam > ${outdir}.stats_summary.txt
	samtools flagstat ${filename}.aln-fixmate.bam > ${outdir}.stats.txt

	sequence="${outdir//stats\/}"
	sequence="${sequence//_L001_R1_001_val_1.fq.gz}"
	
	samtools depth ${filename}.aln-fixmate.bam > ${outdir}.depth.txt
	
	cov=`samtools depth -a ${filename}.aln-fixmate.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'`
	echo $sequence $cov >> stats/CoverageBreadth.txt
	
    dep=`samtools depth -a ${filename}.aln-fixmate.bam | awk '{c++;s+=$3}END{print s/c}'`
    echo $sequence $dep >> stats/MeanReadDepth.txt
    
	nb=`grep ^"SN	raw total sequences:" ${outdir}.stats_summary.txt | cut -d ":" -f2`
	echo $sequence $nb >> stats/NbSeq.txt
	
	qual=`grep ^"SN	average quality:" ${outdir}.stats_summary.txt | cut -d ":" -f2`
	echo $sequence $qual >> stats/Quality.txt

done < list.sam.txt
conda deactivate

conda activate bcftools-1.13
while read line
do
	filename=$line
	filename="${filename//.aln-pe.sam}"
	
	# convert the BAM file into genomic positions. Use mpileup to produce a BCF file that contains all of the locations in the genome. 
	# We use this information to call genotypes and reduce our list of sites to those found to be variant by passing this file into bcftools call.
	bcftools mpileup -Ob -o ${filename}.bcf -f ../snippy/ref/RefSeq_resistance.fasta ${filename}.aln-fixmate.bam     
	bcftools call -vmO z -o ${filename}.vcf.gz ${filename}.bcf
	bcftools call -mO z -o ${filename}.allPos.vcf.gz ${filename}.bcf

	# normalize indels
	bcftools norm -f ../snippy/ref/RefSeq_resistance.fasta ${filename}.vcf.gz -Ob -o ${filename}.norm.vcf.gz

	# filter the data
	bcftools filter -O z -o ${filename}.filtered.vcf.gz -g 12 -i '%QUAL > 50 & MQ > 20 & DP > 20 & (DP4[2]+DP4[3])/sum(DP4) >= 0.9' ${filename}.norm.vcf.gz
    
	# transform as Fasta
	filename2="${filename/mapping_resistance/Fasta}"
	filename3="${filename/mapping_resistance/}"
	filename3="${filename3///}"
	bcftools index $filename.filtered.vcf.gz
	cat ../snippy/ref/RefSeq_resistance.fasta | bcftools consensus -i 'type="snp"' ${filename}.filtered.vcf.gz > ${filename2}.fasta

done < list.sam.txt
conda deactivate

echo 'step 3 complete' > step3.txt

#-g 3:'indel',other
