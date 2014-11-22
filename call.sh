#!/bin/bash

key=$1
lib=$2

#ref="/home/crangen/zding/projects/common/ref/hg38.fa"
ref="/home/crangen/zding/projects/common/ref/old/hg19_GRCh37/hg19.fa"

projectdir="/home/crangen/zding/projects/raindance"
#ID="WTCHG_98544_03"
ID=$key
SM=$ID
#LB="WTCHG_98544"
LB=$lib

samdir="${projectdir}/bams/${ID}"
if [ ! -d  $samdir ]
then
    mkdir ${samdir}
fi

read1="${projectdir}/fastq.qc/${ID}_1.qc.gz"
read2="${projectdir}/fastq.qc/${ID}_2.qc.gz"

sam="${samdir}/${ID}.sam"
bam="${samdir}/${ID}.bam"
sortedbam="${samdir}/${ID}.sorted.bam"

# aligner
bwa mem -R "@RG\tID:$ID\tSM:$SM\tLB:$LB" $ref $read1 $read2 > $sam

samtools fixmate -O bam $sam $bam
samtools sort -T /tmp/${ID}.sorted -o $sortedbam $bam
samtools index $sortedbam

rm $sam $bam

# samtools-exp-rc mpileup -t DP,DV -C50 -Q 30 -d 10000 -u \
#     -l $regions \
#     -b $bamlist \
#     -f $reference \
#     | bcftools-exp-rc call -T $regions -mvO z > $outvcf

# bcftools-exp-rc index $outvcf
