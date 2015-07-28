import sys
import os

import cPickle
import gzip
import argparse
import operator
import numpy as np
import pdb
import re
from Bio import SeqIO

import pysam
import glob

import logging as LG
LG.basicConfig(level=LG.INFO)

params = {
    # wimm
    "sourcedir": "/hts/data6/miseq/agoriely/2015-07-07",
    "fastqdir": "/home/crangen/zding/projects/raindance/mip/fq",
    "bamdir": "/home/crangen/zding/projects/raindance/mip/aln",
    "outdir":"/home/crangen/zding/projects/raindance/mip/scan",
    "design": "/home/crangen/zding/projects/raindance/mip/scan/design.variant",
    # local
    # "sourcedir": "/Users/zd1/volumn/miseq",
    # "fastqdir": "/Users/zd1/cloud/data/raindance/miseq/fq",
    # "bamdir": "/Users/zd1/volumn/wimm/raindance/mip/aln",
    # "outdir":"/Users/zd1/cloud/data/raindance/miseq",
    # "design": "/Users/zd1/cloud/data/raindance/miseq/mip/design.variant",
    "n_randomcode": 5,
    "probelen": 25,
    "ksize": 8
}
    
def run_molid():
    run = Analysis()
    run.count_molid()

        
def samplelist():
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
        
    ofh = open("%s/samples"%params["fastqdir"], "w")
    for tag in tagfq.keys():
        if re.search("Undetermined", tag): 
            continue
        r1fq = tagfq[tag]["R1"]
        r2fq = tagfq[tag]["R2"]
        ofh.write("%s,%s,%s"%(tag, r1fq, r2fq) + "\n")
    ofh.close()
    
def process_fastq():
    samples = []
    read1 = []
    read2 = []
    ofh = open("%s/samples"%params["fastqdir"], "r")
    outlist = open("%s/samples_trim"%params["fastqdir"], "w")
    for line in ofh:
        tag, r1fq, r2fq = line.strip().split(",")
        samples.append(tag)
        read1.append(r1fq)
        read2.append(r2fq)
        outlist.write("%s,%s/%s_1.ps.gz,%s/%s_2.ps.gz\n"%(tag,
                                                       params["fastqdir"], tag,
                                                       params["fastqdir"], tag))
    ofh.close()
    outlist.close()
    
    i = int(os.environ['SGE_TASK_ID']) - 1
    run = Analysis()
    LG.info("scanning run %s, \nread1 %s \nread2 %s "%(samples[i], read1[i], read2[i]))
    run.remove_mid_fastq(read1[i], read2[i], samples[i], params["fastqdir"])


def scan_bam():
    samples = []
    ofh = file("%s/samples"%params["fastqdir"], "r")
    for line in ofh:
        tag = line.strip().split(",")[0]
        samples.append(tag)
    ofh.close()

    design = MIPdesign(params["design"], params["ksize"])
    design.loaddesign()
    design.build_hash()

    i = int(os.environ['SGE_TASK_ID']) - 1
    run = Analysis()
    sam = samples[i]
    bam = "%s/%s.sorted.bam"%(params["bamdir"], sam)
    LG.info("scanning sample %s"%(sam))
    run.classify_bam(sam, bam, design, params["outdir"])
        
class AlleleCounter:
    def __init__(self):
        self.allele = {
            "A":0,
            "C":0,
            "G":0,
            "T":0,
            "N":0}
        
        self.molid = {
            "A":[],
            "C":[],
            "G":[],
            "T":[],
            "N":[]
            }
        
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

    def mapfraction(self, cigar):
        ntotal = 0
        nmap = 0
        for c in cigar:
            if c[0] == 0:
                nmap += c[1]
            ntotal += c[1]
        if ntotal == 0:
            return 0
        else:
            return nmap*1.0/ntotal
        
    def classify_bam(self, sam, bam, design, outdir):
        read1hash = design.read1hash
        read2hash = design.read2hash
        nTotal = 0 # total number of pairs of reads for this run
        nNoMIP = 0 # number of pairs of reads with barcode not found in either read1 or read2
        nNotpaired = 0 # number of pairs of reads with barcode not matched between read1 and read2
    
        samfile = pysam.Samfile("%s"%bam, "rb" )
        ofh = open("%s/%s.count"%(outdir, sam), "w")
        for mip in design.mipvariant.keys():
            mip = "45_FGFR2_exon5_p.C227S_79_2_55"
            variants = design.mipvariant[mip]
            for vt in variants:
                LG.info("mip:%s, variant:%s"%(mip, vt))
                chrm, pos = vt.split(":")
                pos = int(pos)
                ac = AlleleCounter()
                for alignedread in samfile.fetch(chrm, pos-1, pos):
                    nTotal += 1
                    if nTotal % 10000 == 0:
                        LG.info("processed %d reads"%nTotal)
                    ## skip paired reads
                    if not alignedread.is_paired:
                        continue
                    ## skip reads with low fraction of bases mapped
                    if self.mapfraction(alignedread.cigar) < 0.9:
                        continue

                    ## quantify only reads for the current mip
                    seq = alignedread.seq
                    if alignedread.is_read1:
                        mipid, mipanno = design.read1hash.seedkmerMaxOnly(seq[:params["probelen"]], minMatch = 3)
                    else:
                        mipid, mipanno = design.read2hash.seedkmerMaxOnly(seq[:params["probelen"]], minMatch = 3)
                    if mipid != mip:
                        continue
                    
                    mid = alignedread.qname.split("|")[0][1:]
                    
                    # get variant position p
                    pi = [p for p in range(len(alignedread.positions)) if alignedread.positions[p] == pos]
                    if len(pi) == 0:
                        continue
                    nucleotide = alignedread.query[pi[0]-1] # query is 0-based while alignedread is 1-based
                    ac.allele[nucleotide] += 1
                    ac.molid[nucleotide].append(mid)
                    
                row = [sam, mip, vt]
                for b in ["A", "C", "G", "T", "N"]:
                    row.append("%s,%s"%(b, ac.allele[b]))
                    for ndiff in [1,2,3,4]:
                        t = Tag()
                        t.taglist = ac.molid[b]
                        row.append("[%d]%d{%s}"%(ndiff,
                                                 t.uniqtags(ndiff),
                                                 t.histogram(ndiff=2)))
                ofh.write(",".join(row) + "\n")
                del ac

        samfile.close()
        ofh.close()

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

def stringdiff(s1,s2):
    ndiff = 0
    assert(len(s1) == len(s2))
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            ndiff += 1
    return ndiff


class Tag:
    taglist = []
    def __init__(self):
        pass
    
    def uniqtags(self, ndiff = 1):
        tagset = set(self.taglist)
        if ndiff == 1:
            return len(tagset)
        sortedset = sorted(tagset)
        pre = None
        count = 0
        for i in range(len(sortedset)):
            if pre is None:
                count +=1
                pre = sortedset[i]
                continue
            nd = stringdiff(sortedset[i], pre)
            if nd >= ndiff:
                count +=1
                pre = sortedset[i]
        return count

    def histogram(self, ndiff=2, top=10):
        sortedlist = sorted(self.taglist)
        pre = None
        count = 0
        group = {}
        for i in range(len(sortedlist)):
            if pre is None:
                count +=1
                pre = sortedlist[i]
                group[pre] = 1
                continue
            nd = stringdiff(sortedlist[i], pre)
            if not group.has_key(pre):
                group[pre] = 0
            group[pre] += 1
            if nd >= ndiff:
                count +=1
                pre = sortedlist[i]
        sorted_group = sorted(group.items(), key=operator.itemgetter(1), reverse=True)
        topOccur = {}
        for g in sorted_group:
            if not topOccur.has_key(g[1]):
                topOccur[g[1]] = 1
            else:
                topOccur[g[1]] += 1
            if len(topOccur.keys())>= top:
                break
        return sorted(topOccur.items(), reverse=True)
    
def countUniq(inputset, ndiff=1):
    if ndiff == 1:
        return len(inputset)
    sortedset = sorted(inputset)
    pre = None
    count = 0
    for i in range(len(sortedset)):
        if pre is None:
            count +=1
            pre = sortedset[i]
            continue
        nd = stringdiff(sortedset[i], pre)
        if nd >= ndiff:
            count +=1
            pre = sortedset[i]
    return count

def testcountunique():
    x = ["AACCTT",
             "AACCTA",
            "AACCTC",
            "AACCCG",
            "AACCKN"]
    t = Tag()
    t.taglist = x
    print t.uniqtags(ndiff=2)
    print t.histogram(ndiff=2)
    
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
    # testcountunique()
    parser = argparse.ArgumentParser(description='Raindance analysis')
    parser.add_argument('--pairup', action='store_true', default=False)
    parser.add_argument('--fastq', action='store_true', default=False)
    parser.add_argument('--bam', action='store_true', default=False)
    args = parser.parse_args()
    LG.info(args)
    try:
        if args.pairup:
            LG.info("making a sample list")
            samplelist()
        elif args.fastq:
            LG.info("processing fastq")
            process_fastq()
        elif args.bam:
            LG.info("scan bam to quantify alleles")
            scan_bam()
        else:
            print "need to specify --fastq, --bam"
            sys.exit(1)
    except:
        sys.exit(1)
