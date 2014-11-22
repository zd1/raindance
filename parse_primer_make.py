#!/usr/bin/env python3
import os
import sys
from snakemake.utils import read_job_properties

import parse_primer_cfg
params = parse_primer_cfg.param


jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# do something useful with the threads
threads = job_properties[threads]

os.system("qsub -t {threads} {script}".format(threads=threads, script=jobscript))

def getSamples():
    ofh = file(params["samplelist"], "r")
    for line in ofh:
        SAMPLES.append(line.strip())
    ofh.close

SAMPLES = getSamples()

rule all:
    input:
        "/users/mcvean/zd1/volumn/raindance/pileup/{sample}/{sample}.hd5"
    params:
        project = "mcvean.prjc"
        queue =  "short.qc" 
            
rule fastq:
    input:
        "/users/mcvean/zd1/volumn/raindance/fastq/{sample}_1.fastq.gz",
        "/users/mcvean/zd1/volumn/raindance/fastq/{sample}_2.fastq.gz",
    output:
        "/users/mcvean/zd1/volumn/raindance/fastq.qc/{sample}_1.fastq.gz"
        "/users/mcvean/zd1/volumn/raindance/fastq.qc/{sample}_2.fastq.gz"
    shell:
        "python parse_primer.py --fastq {sample}"

rule map:
    input:
        "/users/mcvean/zd1/volumn/raindance/fastq.qc/{sample}_1.fastq.gz"
        "/users/mcvean/zd1/volumn/raindance/fastq.qc/{sample}_2.fastq.gz"
    output:
        "/users/mcvean/zd1/volumn/raindance/bams/${sample}.sorted.bam"
    shell:
        "/users/mcvean/zd1/workspace/rain/map.sh {REF} {read1} {read2} {sample} {library} {outdir}"
        
rule pilup:
    input:
        "/users/mcvean/zd1/volumn/raindance/bams/{sample}.bam"
    output:
        "/users/mcvean/zd1/volumn/raindance/pileup/{sample}/{sample}.hd5"
    shell:
       "python parse_primer.py --pileup {sample}"
