#!/bin/bash -x

# [[file:~/dev/org/pipe.org::*MIP%20analysis][process]]

FILES=$(cat /home/crangen/zding/projects/raindance/mip/scan/*.fastq.gz.molecularid.csv)

for f in "${FILES[@]}"
do
    echo $f
done

# process ends here
