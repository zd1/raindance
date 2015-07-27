#!/bin/bash -x

# [[file:~/dev/org/pipe.org::*MIP%20analysis][alignment]]

module load bwa
module load samtools

ref=/home/crangen/zding/projects/common/ref/old/hg19_GRCh37/hg19.fa
bc=822-1000ng-Probe-1-4000_S26_L001

bwa mem -R "@RG\tID:${bc}\tSM:${bc}\tLB:Test" $ref \
    /hts/data6/miseq/agoriely/2015-05-14reanalysed/${bc}_R1_001.fastq.gz \
    /hts/data6/miseq/agoriely/2015-05-14reanalysed/${bc}_R2_001.fastq.gz > ${bc}.sam

samtools fixmate -O bam .sam ${bc}.bam
samtools sort -T /tmp/${bc}.sorted -o ${bc}.sorted.bam ${bc}.bam
samtools index ${bc}.sorted.bam

# alignment ends here
