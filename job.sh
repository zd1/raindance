#!/bin/bash
 
#$ -N rain
#$ -q short.qc
#$ -P mcvean.prjc
#$ -V
#$ -wd /users/mcvean/zd1/volumn/raindance/log/
#$ -e /users/mcvean/zd1/volumn/raindance/log/$JOB_NAME-$JOB_ID-$TASK_ID.e
#$ -o /users/mcvean/zd1/volumn/raindance/log/$JOB_NAME-$JOB_ID-$TASK_ID.o
#$ -t 3-288

python /users/mcvean/zd1/workspace/raindance/parse_primer_make.py





