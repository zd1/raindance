#!/bin/bash

${HOME}/opt/vpythons/3.3/bin/snakemake --directory ${HOME}/volumn/raindance/log/ \
       --snakefile ${HOME}/workspace/raindance/parse_primer_make.py \
       --forceall \
       --cluster "qsub -P mcvean.prjc -q short.qc" \
       -j 2 -p \
       all
