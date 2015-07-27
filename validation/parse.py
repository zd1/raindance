
import sys
import os

import cPickle
import gzip
import argparse
import numpy as np
import pdb
import re
from Bio import SeqIO

import pysam
import glob

import logging as LG
LG.basicConfig(level=LG.INFO)

params = {
    "sourcedir": "/hts/data6/miseq/agoriely/2015-07-07",
    "fastqdir": "/home/crangen/zding/projects/raindance/mip/fq",
    "bamdir": "/home/crangen/zding/projects/raindance/mip/bam",
    "outdir":"/home/crangen/zding/projects/raindance/mip/scan",
    "design": "/home/crangen/zding/projects/raindance/mip/scan/design.variant",

    # "sourcedir": "/Users/zd1/volumn/miseq",
    # "fastqdir": "/Users/zd1/cloud/data/raindance/miseq/fq",
    # "bamdir": "/Users/zd1/cloud/data/raindance/miseq/bam",
    # "outdir":"/Users/zd1/cloud/data/raindance/miseq",
    # "design": "/Users/zd1/cloud/data/raindance/miseq/mip/design.variant",
    "n_randomcode": 5,
    "probelen": 30,
    "ksize": 8
}
    
def run_molid():
    run = Analysis()
    run.count_molid()

def run_bam():
    design = MIPdesign(params["design"], params["ksize"])
    design.loaddesign()
    design.build_hash()
    
    run = Analysis()
    bams = glob.glob("%s/*.sorted.bam"%params["sourcedir"])
    for bam in bams:
        m = re.search(".*\/(.*)_(S.+)_(L.+).sorted.bam", bam)
        tag,sam,ln = m.groups()
        LG.info("scanning run %s, sample %s"%(tag, sam ))
        run.classify_bam(bam, tag, design, params["outdir"])
            
def process_fastq():
    design = MIPdesign(params["design"], params["ksize"])
    design.loaddesign()
    design.build_hash()
    
    run = Analysis()
    fastqs = glob.glob("%s/*fastq.gz"%params["sourcedir"])

    tagfq = {}
    for fq in fastqs:
        m = re.search(".*\/(.*)_(S.+)_(L.+)_(R.+)_001.fastq.gz",fq)
        tag,sam,ln,rd = m.groups()
        if not tagfq.has_key(tag):
            tagfq[tag] = {}
            tagfq[tag]["Lane"] = []
            tagfq[tag]["Sam"] = []
        tagfq[tag][rd]=fq
        tagfq[tag]["Lane"].append(ln)
        tagfq[tag]["Sam"].append(sam)
    
    for tag in tagfq.keys():
        if re.search("Undetermined", tag): 
            continue
        r1fq = tagfq[tag]["R1"]
        r2fq = tagfq[tag]["R2"]
        LG.info("scanning run %s, \nread1 %s \nread2 %s "%(tag, r1fq, r2fq ))
        run.remove_mid_fastq(r1fq, r2fq, tag, params["fastqdir"])

class BCcounter:
    def __init__(self):
        self.nTotal = 0
        self.r1basecount = {
            "A":0,
            "C":0,
            "G":0,
            "T":0,
            "N":0,
            }
        self.r2basecount = {
            "A":0,
            "C":0,
            "G":0,
            "T":0,
            "N":0,
            }

        self.randoms = set([])

        
class MIPdesign:
    '''note we expect the design file in the following format, 
    48_FGFR2_exon4_p.Q174P_82_2_57,AAGGTGAGGACTTTCTGAATCTA,GTCGGAGGAGACGTAGA,123279473,123279624,-,chr10:123279605
    '''
    mipid = []
    read1prob = []
    read2prob = []
    read1prob_dict = {} # seq -> id
    read2prob_dict = {} # seq -> id
    read1hash = {} # end 
    read2hash = {} # start

    # 7th column is for the variant position
    # [mipid][variant pos]
    mipvariant = {}
    mippos = {}
    ksize = None
    def __init__(self, designtemplate, ksize):
        self.designtemplate = designtemplate
        self.ksize = ksize
        pass
    def loaddesign(self):
        head=True
        ofh = file(self.designtemplate, 'r')
        for line in ofh:
            if head:
                head = False
                continue
            d=line.strip().split(",")
            LG.debug("porbes:%s"%str(d))
            pid = d[0] # ID column
            variantpos = d[6] # variant position chrx:pos
            mipstart = int(d[3]) # mip start position
            mipend = int(d[4]) # mip end position
            self.mipid.append(pid)
            if not self.mipvariant.has_key(pid):
                self.mipvariant[pid]=[variantpos]
            else:
                self.mipvariant[pid].append(variantpos)
            if not self.mippos.has_key(pid):
                self.mippos[pid] = [mipstart, mipend]
            read1 = revecomp("".join(d[1].split()).upper())
            read2 = "".join(d[2].split()).upper()
            self.read1prob.append(read1)
            self.read2prob.append(read2)
            self.read1prob_dict[read1] = pid
            self.read2prob_dict[read2] = pid
        ofh.close()
        assert len(self.mipid) == len(self.read1prob) and len(self.mipid) == len(self.read2prob)
        LG.info("loaded %d probes"%(len(self.mipid)))

    def build_hash(self):
        self.read1hash = SeqHash(params["ksize"], self.read1prob_dict)
        self.read2hash = SeqHash(params["ksize"], self.read2prob_dict)

class Analysis:
    bams = []
    barcodes_seq_dict = {}
    barcodes_dict = {}
    bchash = None

    def __init__(self):
        pass

    def _barcodebasematrix(self):
        return np.zeros((5,params["barcodearea"])) # matrix for base frq for the first 6 bases

    def _baseidx(self, base):
        if base == 'A':
            return 0
        elif base == 'T':
            return 1
        elif base == 'G':
            return 2
        elif base == 'C':
            return 3
        else:
            return 4

    def _baseidx_reverse(self, idx):
        if idx == 0:
            return 'A'
        elif idx == 1:
            return 'T'
        elif idx == 2:
            return 'G'
        elif idx == 3:
            return 'C'
        elif idx == 4:
            return 'N'

    def classify_bam(self, bam, tag, design, outdir):
        read1hash = design.read1hash
        read2hash = design.read2hash
        nTotal = 0 # total number of pairs of reads for this run
        nNoMIP = 0 # number of pairs of reads with barcode not found in either read1 or read2
        nNotpaired = 0 # number of pairs of reads with barcode not matched between read1 and read2
        bcounter = {}
        bcounter[tag] = {}
    
        samfile = pysam.Samfile("%s"%bam, "rb" )

        for mip in design.mipvariant.keys():
            if not bcounter[tag].has_key(mip):
                bcounter[tag][mip] = {}

            variants = design.mipvariant[mip]
            for vt in variants:
                chrm, pos = vt.split(":")
                pos = int(pos)
                if not bcounter[tag][mip].has_key(pos):
                    bcounter[tag][mip][pos] = BCcounter()
                posdata = {}
                rname_mid = {}
                for alignedread in samfile.fetch(chrm, pos-1, pos):
                    nTotal += 1
                    if nTotal % 10000 == 0:
                        LG.info("processed %d reads"%nTotal)

                    seq = alignedread.seq
                    if alignedread.is_read1:
                        mipid, anno = read1hash.seedkmerMaxOnly(seq[params["n_randomcode"]:params["probelen"]], 3)
                    else:
                        mipid, anno = read2hash.seedkmerMaxOnly(seq[params["n_randomcode"]:params["probelen"]], 3)
                    if mip != mipid:
                        nNoMIP +=1
                        continue

                    if not posdata.has_key(alignedread.qname):
                        # bases by [r1,r2]
                        posdata[alignedread.qname] = np.zeros((5,2), 'i8')
                    if not rname_mid.has_key(alignedread.qname):
                        rname_mid[alignedread.qname] = {"r1":"", "r2":""}

                    if alignedread.is_read1:
                        rname_mid[alignedread.qname]["r1"] = seq[:params["n_randomcode"]]
                        ri = 0
                    else:
                        rname_mid[alignedread.qname]["r2"] = seq[:params["n_randomcode"]]
                        ri = 1
                    
                    # get variant position p
                    pi = [p for p in range(len(alignedread.positions)) if alignedread.positions[p] == pos]
                    if len(pi) == 0:
                        continue

                    bi = basei(alignedread.query[pi[0]])
                    posdata[alignedread.qname][bi,ri] += 1

                for readpair in posdata.keys():
                    _mid = rname_mid[readpair]["r1"] + rname_mid[readpair]["r2"]
                    if len(_mid) < 10:
                        nNotpaired += 1 
                    bcounter[tag][mip][pos].randoms.add(_mid)
                    bcounter[tag][mip][pos].nTotal = np.sum(posdata[readpair])
                    for b in ["A", "C", "G", "T", "N"]:
                        bcounter[tag][mip][pos].r1basecount[b] = posdata[readpair][basei(b), 0]
                        bcounter[tag][mip][pos].r2basecount[b] = posdata[readpair][basei(b), 1]
                break
        samfile.close()
        
        ofh = open("%s/%s.count"%(outdir, tag), "w")        
        for p in bcounter[tag][mip].keys(): # for each variant position
            row = [tag, mip, nTotal, nNoMIP, nNotpaired,
                   bcounter[tag][mip][p].nTotal, len(bcounter[tag][mip][p].randoms)]
            pdb.set_trace()
            for r in [0,1]:
                if r == 0:
                    ct = bcounter[tag][mip][pos].r1basecount
                else:
                    ct = bcounter[tag][mip][pos].r2basecount
                for b in ["A", "C", "G", "T", "N"]:
                    row.append("%s:%d"%(b, ct[b]))
            ofh.write(",".join(row) + "\n")

        ofh.close()
                    
        pass

    def remove_mid_fastq(self, r1fq, r2fq, tag, outdir):
        '''remote molcular IDs from reads in fastq files.
        the id is consisted of the first five bases of read1 plus
        the first five bases of read2. we add that to the read name.
        Note that we need to ensure that read1 and read2 have a same
        name.
        '''
        fastq_read1 = r1fq
        fastq_read2 = r2fq
        LG.info("read1 file %s "%fastq_read1)
        LG.info("read2 file %s "%fastq_read2)
        
        fastq_out_read1 = outdir + "/" + tag  + "_1.ps.gz"
        fastq_out_read2 = outdir + "/" + tag  + "_2.ps.gz"

        LG.info("read1 output file %s "%fastq_out_read1)
        LG.info("read2 output file %s "%fastq_out_read2)
        
        LG.info("%s"%tag)
        LG.info(" read1:%s"%fastq_read1)
        LG.info(" read2:%s"%fastq_read2)
        
        handle1 = gzip.open(fastq_read1, "r")
        handle2 = gzip.open(fastq_read2, "r")
        output_handle1 = gzip.open(fastq_out_read1, "w")
        output_handle2 = gzip.open(fastq_out_read2, "w")
            
        record1_iterator = SeqIO.parse(handle1, "fastq")
        record2_iterator = SeqIO.parse(handle2, "fastq")
        
        start = params["n_randomcode"]

        c = 0
        while True:
            c+=1
            if (c%10000 == 0):
                LG.info("scanned %d reads"%c)
            try:
                record1 = next(record1_iterator, None)
                record2 = next(record2_iterator, None)
            except:
                LG.error("error when scanning, skip read")
                continue
            if not record1 or not record2:
                LG.info("skip read read1:%s read2:%s"%(record1, record2))
                break

            mid = str(record1.seq[:start]) + str(record2.seq[:start])
            record1 = record1[start:]
            record2 = record2[start:]
            record1.id = mid + "|" + record1.id
            record2.id = mid + "|" + record2.id
            record1.description = ""
            record2.description = ""
            output_handle1.write(record1.format("fastq"))
            output_handle2.write(record2.format("fastq"))

        handle1.close()
        handle2.close()
        output_handle1.close()
        output_handle2.close()

        LG.info("Completed scanning for %s"%tag)

    
    def classify_fastq(self, r1fq, r2fq, tag, design, outdir):
        
        read1hash = design.read1hash
        read2hash = design.read2hash
        nTotal = 0 # total number of pairs of reads for this run
        nNoMIP = 0 # number of pairs of reads with barcode not found in either read1 or read2
        nNotpaired = 0 # number of pairs of reads with barcode not matched between read1 and read2
        bcounter = {}
        # count base frequencies
        molidmx = np.zeros((5,2*params["n_randomcode"]))
        LG.info(" tag:%s"%tag)
        with gzip.open(r1fq, 'r') as hdl1, gzip.open(r2fq, 'r') as hdl2, open("%s/%s.molid"%(outdir, tag), "w") as ofhid:
            for read1, read2 in zip(SeqIO.parse(hdl1, "fastq"), SeqIO.parse(hdl2, "fastq")):
                nTotal += 1
                if nTotal % 10000 == 0:
                    LG.info("processed %d reads"%nTotal)
                mipid1, anno1 = read1hash.seedkmerMaxOnly(str(read1.seq)[params["n_randomcode"]:params["probelen"]], 1)
                mipid2, anno2 = read2hash.seedkmerMaxOnly(str(read2.seq)[params["n_randomcode"]:params["probelen"]], 1)
                if mipid1 is None or mipid2 is None:
                    nNoMIP += 1
                    continue
                if mipid1 != mipid2:
                    nNotpaired += 1
                    continue
                if not design.mipvariant.has_key(mipid1):
                    LG.error("found mip ID not in design %s"%mipid1)
                    continue
                

                # tag -> mips -> molids
                if not bcounter.has_key(tag):
                    bcounter[tag] = {}
                if not bcounter[tag].has_key(mipid1):
                    bcounter[tag][mipid1] = {"read1": {}, "read2": {}}
                    
                # molecular ID. combine r1 and r2 as they are from a same molecule
                molid = str(read1.seq)[:params["n_randomcode"]] + str(read2.seq)[:params["n_randomcode"]]

                mipstart, mipend = design.mippos[mipid1]
                offsetextend = 4
                for vpos in design.mipvariant[mipid1]:
                    assert mipstart < vpos and mipend > vpos, "error in variant position %d"%vpos
                    offset1 = vpos - mipstart + params["n_randomcode"]
                    if offset1 >= 0 and offset1 < 100:
                        vbase_read1 = str(read1.seq)[(offset1-offsetextend):(offset1+offsetextend+1)]
                        if not bcounter[tag][mipid1]["read1"].has_key(vpos):
                            bcounter[tag][mipid1]["read1"][vpos]={}
                        if not bcounter[tag][mipid1]["read1"][vpos].has_key(vbase_read1):
                            bcounter[tag][mipid1]["read1"][vpos][vbase_read1] = BCcounter()
                        bcounter[tag][mipid1]["read1"][vpos][vbase_read1].nTotal +=1
                        bcounter[tag][mipid1]["read1"][vpos][vbase_read1].randoms.add(molid)
                        
                    offset2 = mipend - vpos + params["n_randomcode"]
                    vbase_read2 = str(read2.seq)[(offset2-offsetextend):(offset2+offsetextend+1)]
                    if offset2 >= 0 and offset2 < 100:
                        if not bcounter[tag][mipid1]["read2"].has_key(vpos):
                            bcounter[tag][mipid1]["read2"][vpos]={}
                        if not bcounter[tag][mipid1]["read2"][vpos].has_key(vbase_read2):
                            bcounter[tag][mipid1]["read2"][vpos][vbase_read2] = BCcounter()
                        bcounter[tag][mipid1]["read2"][vpos][vbase_read2].nTotal +=1
                        bcounter[tag][mipid1]["read2"][vpos][vbase_read2].randoms.add(molid)
                # if mipid1 in {"92_PTEN_exon8_p.T319P_49_2_46", "431_FGFR1_exon12_p.K477N_600_2_21",
                #               "131_HRAS_exon3_p.Q43H_166853_2_63", "222_TP53_exon2_p.G67V_0103_3_30"}:
                #     ofhid.write("%s,%s\n"%(mipid1, molid))

                
                ## read1.seq (find the variant A/G)
                ## read2.seq (find the mutation)
                ## print "%s,%s,%s,%s"%(tag, molid, read1.seq[52:55], read2.seq[24:30])

                for bi in range(len(molid)):
                    molidmx[self._baseidx(molid[bi]),bi] += 1
                
        with open("%s/%s.count"%(outdir, tag), "w") as ofh:
            for mip in bcounter[tag].keys(): # for each sample
                for rd in ["read1", "read2"]: # for each read
                    if bcounter[tag][mip].has_key(rd):
                        for p in bcounter[tag][mip][rd].keys(): # for each variant position
                            record = bcounter[tag][mip][rd][p]
                            for nts in record.keys(): # write out pos, number of reads, number of molIDs
                                row = [tag, mip, nTotal, nNoMIP, nNotpaired]
                                row.append("rd:%s|nt:%s|p:%s|t:%d|r:%d"%(rd, nts, p,record[nts].nTotal, len(record[nts].randoms)))
                                row = map(str, row)
                                ofh.write(",".join(row) + "\n")
                
        # convert to frequency
        molidmx = molidmx/molidmx.sum(axis=0)
        with open("%s/%s.basemx"%(outdir, tag), "w") as ofh:
            for i in range(5):
                row = [tag]
                row.extend(molidmx[i,:])
                ofh.write(",".join(map(str, row))+ "\n")

    def compare_str(self, str1, str2):
        ndiff  = 0
        assert len(str1) == len(str2), "the length of str1 doesn't equal to that of str2"
        for i in range(len(str1)):
            if str1[i] != str2[i]:
                ndiff += 1
        return ndiff
                
    def count_molid(self):
        '''we want to compare molids and count the number of basepairs
        that differ between two neighbouring IDs. 
        Here we assume the list is sorted.
        the first column of the input file should be MIP ID, and the
        second column should be molid - 10 basepair nuleotides
        '''
        # [MIP ID][Mol ID] = number
        count = {}
        previous_molid = ""
        with file("/home/crangen/zding/projects/raindance/mip/scan/Tes3D2.molid.sorted", "r") as ofh:
            for line in ofh:
                item = line.strip().split(",")
                mipid, molid = item
                if not count.has_key(mipid):
                    count[mipid] = {}
                if not count[mipid].has_key(molid):
                    count[mipid][molid] = 0
                if previous_molid == "":
                    count[mipid][molid] = 0
                elif previous_molid == molid:
                    pass
                else:
                    ## we want to compare the current molid against the previous one
                    ndiff = self.compare_str(molid, previous_molid)
                    count[mipid][molid] = ndiff
                previous_molid = molid
                
        with file("/home/crangen/zding/projects/raindance/mip/scan/Tes3D2.molid.diff", "w") as ofh:
            for mip in count.keys():
                for mid in count[mip].keys():
                    row = [mip, mid, str(count[mip][mid])]
                    ofh.write(",".join(row) + "\n")
            pass


class SeqHash:
    '''build a kmer hash'''
    class Kmer:
        def __init__(self, size, pos, seq, tag):
            self.size = size
            self.pos = pos
            self.seq = seq
            self.tag = tag

    def __init__(self, kmersize, sequence_dict):
        self.kmersize = kmersize
        self._buildHash(sequence_dict)

    def _buildHash(self, sequence_dict):
        '''sequence_dict contains seq -> kmer object'''
        resulthash = {}
        for seq in sequence_dict.keys():
            seqlen = len(seq)
            tag = sequence_dict[seq]
            if seqlen < self.kmersize:
                LG.error("can't use kmer that"
                         " is larger then the seuqence %s (kmersize%d)"%(seq, self.kmersize))
                sys.exit(1)
            for i in range(seqlen - self.kmersize):
                start = i
                end = i + self.kmersize
                kmer = self.Kmer(self.kmersize, i, seq[start:end], tag)
                if not resulthash.has_key(kmer.seq):
                    resulthash[kmer.seq] = []
                resulthash[kmer.seq].append(kmer)
        self.hashtable = resulthash


    def read2kmer(self, read, ksize, tag="read"):
        assert len(read) > ksize
        kmers = []
        for i in range(len(read) - ksize):
            start = i
            end = i + ksize
            seq = read[start:end]
            kmer = self.Kmer(ksize, i, seq, tag)
            kmers.append(kmer)
        return kmers

    def seedkmer(self, read):
        kmers = self.read2kmer(read, self.kmersize)
        ids={}
        for kmer in kmers: # read kmers
            if self.hashtable.has_key(kmer.seq): # hashed kmers
                for ksam in self.hashtable[kmer.seq]:
                    if not ids.has_key(ksam.tag):
                        ids[ksam.tag] = []
                    ids[ksam.tag].append([kmer.pos, ksam.pos])
        return ids

    def seedkmerMaxOnly(self, read, minMatch=3):
        '''consider the one with maximum amount of kmer match the match for
        the read'''
        allmatches = self.seedkmer(read)
        n_most = 0
        best=None
        sam = None
        if len(allmatches.keys()) > 0:
            for key in allmatches.keys():
                if len(allmatches[key]) > n_most:
                    n_most = len(allmatches[key])
                    best = allmatches[key]
                    sam = key
        if n_most > minMatch:
            return sam, best
        else:
            return None, None

def cdmz(myobject, filename, protocol = -1):
        """Saves a compressed object to disk
        """
        myfile = gzip.open(filename, 'wb')
        myfile.write(cPickle.dumps(myobject, protocol))
        myfile.close()

def clz(filename):
        """Loads a compressed object from disk
        """
        myfile = gzip.GzipFile(filename, 'rb')
        mybuffer = ""
        while True:
                data = myfile.read()
                if data == "":
                        break
                mybuffer += data
        myobject = cPickle.loads(mybuffer)
        myfile.close()
        return myobject


def revecomp(input):

    alphabet = {
        'A' : 'T',
        'G' : 'C',
        'T' : 'A',
        'C' : 'G',
        'N':"N"}
    output = ""
    for i in range(len(input))[::-1]:
        output += alphabet[input[i]]
    return output

def basei(basestr):
    alphabet = {
        'A' : 0,
        'C' : 1,
        'G' : 2,
        'T' : 3,
        'N': 4}
    return alphabet[basestr]
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Raindance analysis')
    parser.add_argument('--fastq', action='store_true', default=False)
    args = parser.parse_args()
    LG.info(args)
    try:
        if args.fastq:
            LG.info("processing fastq")
            process_fastq()
        else:
            print "need to specify --fastq, --bam"
            sys.exit(1)
    except:
        sys.exit(1)
