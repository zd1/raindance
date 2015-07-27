#!/bin/bash -x

# [[file:~/dev/org/pipe.org::*MIP%20analysis][parse_miseq_job]]

#$ -N miseq
#$ -V
#$ -wd /home/crangen/zding/projects/raindance/mip/scan
#$ -e /home/crangen/zding/projects/raindance/mip/scan/$JOB_NAME.$JOB_ID-$TASK_ID.e
#$ -o /home/crangen/zding/projects/raindance/mip/scan/$JOB_NAME.$JOB_ID-$TASK_ID.o

module load python
/usr/bin/time python parse.py

# parse_miseq_job ends here
