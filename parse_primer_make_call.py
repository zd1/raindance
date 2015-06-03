
import os
import sys
sys.path.append("/users/mcvean/zd1/workspace/raindance/")

import argparse
import csv

import parse_primer_cfg
params = parse_primer_cfg.param

def load_probe_annotation():
    '''read amplicon annotation data, load this file to mem
    probes are identified by their names'''
    seen = set()
    with open(params['table'], 'rt') as csvfile:
        head=True
        pbreader = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in pbreader:
            if head:
                head=False
                continue
            tag = row[0] # we expect this to be unique name, here we use the Amplicon ID
            gene = row[1]
            region = row[11].split() # this is the target region
            chrm = region[0]
            start = int(region[1])
            end = int(region[2])
            ampstart = int(row[22].replace(",",""))
            ampend = int(row[23].replace(",",""))
            ampkey = "%s:%s:%s:%s-%s"%(tag, gene, chrm, start, end)
            if ampkey not in seen:
                seen.add(ampkey)
            else:
                print("Duplicated identifier:%s"%tag)
                sys.exit()
    return [s for s in seen]

def sum_calls_by_region():
    
    pass

probes = load_probe_annotation()

code = "/users/mcvean/zd1/volumn/code/raindance"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Raindance job runner')
    args = parser.parse_args()
    allindics = range(len(probes))
    
    # jobidx
    # index = int(os.environ['SGE_TASK_ID']) - 1
    # batch = 10
    # start = index * batch
    # end = start + batch
    start = 0
    end = 500
    for i in range(start,end):
        probe = probes[allindics[i]]
        print "Analysing probe %s"%probe
        if os.path.exists("/users/mcvean/zd1/volumn/raindance/sum/%s/%s.A.csv"%(probe, probe)):
            print "probe already done"
            continue
                
        # calling for probes
        cmd = "python %s/parse_primer.py"%code
        #cmd += " --agg "
        cmd += " --call "
        cmd += " --amp " + probe
        print cmd
        os.system(cmd)
