#!/bin/bash
 
#$ -N pileup
#$ -q batchq
#$ -V
#$ -e /home/crangen/zding/projects/raindance/log.pileup/$JOB_NAME-$JOB_ID-$TASK_ID.e
#$ -o /home/crangen/zding/projects/raindance/log.pileup/$JOB_NAME-$JOB_ID-$TASK_ID.o

module load python/2.7.5
module load bwa

python /home/crangen/zding/workspace/rain/parse_primer.py 





