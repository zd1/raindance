#!/bin/bash -x

# [[file:~/dev/org/rain.org::*Process][run_job_tangle_local]]

#$ -N fq
#$ -j y
#$ -V
#$ -e /home/crangen/zding/projects/raindance/mip/
#$ -o /home/crangen/zding/projects/raindance/mip/
#$ -t 1-8

module load python
module load bwa
module load samtools

##########  split reads ##########

# python /home/crangen/zding/projects/raindance/mip/code/parse.py --fastq

########## alignment  ##########

# ref=/home/crangen/zding/projects/common/ref/old/hg19_GRCh37/hg19.fa
# data=/hts/data6/miseq/agoriely/2015-07-07
# fqdir=/home/crangen/zding/projects/raindance/mip/fq
# outdir=/home/crangen/zding/projects/raindance/mip/aln
# samples=($(cat ${fqdir}/samples_trim))

# lib=$data

# cd $outdir

# index=$(($SGE_TASK_ID-1))

# sam=$(echo ${samples[$index]} | cut -d, -f1)
# read1=$(echo ${samples[$index]} | cut -d, -f2)
# read2=$(echo ${samples[$index]} | cut -d, -f3)

# echo "processing $sam $read1 $read2"

# bwa mem -R "@RG\tID:${sam}\tSM:${sam}\tLB:${lib}" $ref $read1 $read2 > ${sam}.sam

# samtools fixmate -O bam ${sam}.sam ${sam}.bam
# samtools sort -T /tmp/${sam}.sorted -o ${sam}.sorted.bam ${sam}.bam
# samtools index ${sam}.sorted.bam

########## quantify allele  ##########
python /home/crangen/zding/projects/raindance/mip/code/parse.py --bam

# run_job_tangle_local ends here
