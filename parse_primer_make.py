
import os
import sys
sys.path.append("/users/mcvean/zd1/workspace/raindance/")

import argparse

import parse_primer_cfg
params = parse_primer_cfg.param

SAMPLES = [] 
ofh = open(params["samplelist"], "r")
for line in ofh:
    SAMPLES.append(line.strip())
ofh.close

code = "/users/mcvean/zd1/volumn/code/raindance"

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Raindance job runner')
    parser.add_argument('--indexlist', required=False, default= None)
    args = parser.parse_args()
    allindics = []
    if args.indexlist is not None:
        ofh = open(args.indexlist)
        for line in ofh:
            d = line.strip()
            allindics.append(d)
        ofh.close()
        print allindics
    else:
        allindics = range(len(SAMPLES))

    # jobidx
    index = int(os.environ['SGE_TASK_ID']) - 1
    sample = SAMPLES[allindics[index]]
    print "Analysing sample %s"%sample
    
    wtchg, lib, libsam = sample.split("_")

    # fastq
    if not os.path.exists("%s/%s_1.qc.gz"%(params["fastq_outdir"],sample)) or not os.path.exists("%s/%s_2.qc.gz"%(params["fastq_outdir"],sample)):
            
        cmd = "python %s/parse_primer.py"%code
        cmd += " --fastq "
        cmd += " --sample " + sample
        print cmd
        os.system(cmd)

    # map
    if not os.path.exists("%s/%s/%s.sorted.bam"%(params["bamdir"], sample, sample)):
        outdir = "%s/%s"%(params["bamdir"], sample)
        os.mkdir(outdir)
                
        cmd = "/bin/bash %s/map.sh"%code
        cmd += " " + params["ref"]
        cmd += " " + params["fastq_outdir"]+"/"+sample+"_1.qc.gz"
        cmd += " " + params["fastq_outdir"]+"/"+sample+"_2.qc.gz"
        cmd += " " + sample # sample
        cmd += " " + wtchg + "_" + libsam # library
        cmd += " " + outdir
        print cmd
        os.system(cmd)

    # pileup
    cmd = "python %s/parse_primer.py"%code
    cmd += " --pileup "
    cmd += " --sample " + sample
    #cmd += " --amplist /users/mcvean/zd1/workspace/raindance/amp_shortlist"
    print cmd
    os.system(cmd)
