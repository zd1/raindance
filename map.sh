#!/bin/bash

ref=$1
read1=$2
read2=$3
sample=$4
library=$5
outdir=$6

ID=$sample
SM=$sample
LB=$library

sam="${outdir}/${sample}.sam"
bam="${outdir}/${sample}.bam"
sortedbam="${outdir}/${sample}.sorted.bam"

# aligner
bwa mem -R "@RG\tID:$ID\tSM:$SM\tLB:$LB" $ref $read1 $read2 > $sam

samtools fixmate -O bam $sam $bam
samtools sort -T /tmp/${ID}.sorted -o $sortedbam $bam
samtools index $sortedbam

rm $sam $bam
