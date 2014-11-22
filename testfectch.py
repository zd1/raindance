import sys
import os
if "MYPYTHONPATH" in os.environ:
    sys.path.append(os.environ['MYPYTHONPATH'])
else:
    sys.path.append("/Users/zd1/cloud/myworkspace/tools/src")

import numpy as np
import logging as LG
LG.basicConfig(level=LG.DEBUG)

import pysam
import pdb
import re

params = {"baseQ_cutoff": 20,
          "mapQ_cutoff": 20 # 0-255, -10log10{mapping position is wrong}
          }

def letterpos(seq):
    basepos = {"A":0,"T":1,"G":2,"C":3, "N":4}
    idx = []
    for b in seq:
        k = b.upper()
        if basepos.has_key(k):
            idx.append(basepos[k])
        else:
            idx.append(4)
    return idx

def fastqQuality(ascii_str):
    return [ord(i)-33 for i in ascii_str]

def cigarstr2seq(cigar):
    cigar = re.findall(r'([0-9]+)([MIDNSHPX=])', seq)
    cigar_read = ""
    for i in cigar:
        count, sym = i
        count = int(count)
        for j in range(count):
            cigar_read += sym
    return cigar_read

def cigartuple2seq(tuples):
    seq = ""
    for t in tuples:
        op, c = t
        seq += str(op)*c
    return seq

if __name__ == '__main__':
    
    bam = "/Users/zd1/mount/wimmhts/raindance/bams/WTCHG_98544_01/WTCHG_98544_01.sorted.bam"
#    bam = "test.1.bam"
    samfile = pysam.Samfile("%s"%bam, "rb")
    samid = "test"
    # K-ras
    chrm= "chr1"
    p0 = 17185
    p1 = 17185
    LG.info('%s %s:%d-%d' %(samid, chrm, p0, p1))
    MM = np.zeros((p1-p0,5))
    QQ = np.zeros((p1-p0,2)) # store quality scores, first dimension for read1.
    print MM.shape
    nTotal= 0
    
    for alignedread in samfile.fetch(chrm, p0-1, p1+1):
        # alignedread.qual 
        nTotal += 1
        rleft = None
        rright = None
        seq = None
        print alignedread.cigar
        print alignedread.cigarstring
        print cigartuple2seq(alignedread.cigar)
        
        #  ----t-----
        #------r--------
        if alignedread.pos < p0 and alignedread.pos + alignedread.rlen >= p1:
            rleft = 0
            rright = p1-p0
            seq = alignedread.seq[(p0-alignedread.pos):(p1 - alignedread.pos)]
        #  ----t-----
        #------r---
        elif alignedread.pos < p0 and alignedread.pos + alignedread.rlen < p1:
            rleft = 0
            rright = alignedread.pos + alignedread.rlen - p0
            seq = alignedread.seq[(p0-alignedread.pos):]
        #  ----t-----
        #   ---r---
        elif alignedread.pos >= p0 and alignedread.pos + alignedread.rlen < p1:
            rleft = alignedread.pos - p0
            rright =  alignedread.pos + alignedread.rlen - p0
            seq = alignedread.seq
        #  ----t-----
        #   ---r------
        elif alignedread.pos >= p0 and  alignedread.pos + alignedread.rlen >= p1:
            rleft = alignedread.pos - p0
            rright = p1 - p0
            seq = alignedread.seq[:(p1-alignedread.pos)]
        else:
            continue

        #print "%d - %d (alignreadpos %d, readlength %d):  %s"%(rleft,rright, alignedread.pos, alignedread.rlen, seq)
        idx = letterpos(seq)
        print idx
        for i in range(len(idx)):
            MM[rleft+i,idx[i]] += 1
        if nTotal > 1000:
            break
        
    print MM
    print nTotal
    
    samfile.close()
