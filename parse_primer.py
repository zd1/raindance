import sys
import os
if "MYPYTHONPATH" in os.environ:
    sys.path.append(os.environ['MYPYTHONPATH'])
else:
    sys.path.append("/Users/zd1/cloud/myworkspace/tools/src")

import numpy as np
import h5py
import pdb
import csv
from myio import serial
import gzip
import re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
import parse_primer_cfg

import matplotlib
matplotlib.use('PDF')
import pylab
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pysam
import subprocess
import argparse
import hashlib

import logging as LG
LG.basicConfig(level=LG.DEBUG)

import cProfile, pstats, StringIO

params = parse_primer_cfg.param

def plotAllFastq():
    runner = Analysis()
    runner.load_seq_index()
    runner.load_probe_annotation()
    runner.globalsummary()

def run_pileup():
    runner = Pileup()
    runner.countpileup()
    
def run():

    # pr = cProfile.Profile()
    # pr.enable()
    # # ... do something ...
    
    # index = int(os.environ['SGE_TASK_ID']) - 1
    # runner = Analysis()
    # runner.load_seq_index()
    # runner.load_probe_annotation()
    # runner.build_hash()
    # runner.test_hash()
    # runner.scan_fastq(index)
    # runner.summary_plots(index)
    # runner.run_aligner(index)
    # run_pileup(index)
    
    # pr.disable()
    # s = StringIO.StringIO()
    # sortby = 'cumulative'
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print s.getvalue()

    
    pass
    
class Analysis:
    '''run analysis '''
    # pairwise alignment matching scores
    matchingscore = 2
    mismatchscore = -1
    gapopenpanelty = -1
    gapextendpanelty = -0.1

    # variable to store all probe information. identified by forward
    # sequence. fwseq -> Probe object. 
    probeset = {}

    # fastq list, indexed by libname + barcode
    fastqset = {}
    fastqtag = [] # keys only in same order as input
    
    # SeqHash object, built for primer sequences
    primerHash = None
    
    def __init__(self):
        pass

    def load_seq_index(self):
        '''load fastq index into fastqset dict
        indexed by library ID and barcode, e.g.
        WTCHG_107218_01
        dict in structure:
        {"WTCHG_107218_01":{"1": read1; "2":read2}}'''
        
        with open(params['fastqindex'], 'rt') as csvfile:
            head=True
            seqreader = csv.reader(csvfile, delimiter=',', quotechar='"')
            for row in seqreader:
                # /hts/data3/agoriely/data_Mar2014/orig_data/WTCHG_107218_01_2.fastq.gz
                fastq =  row[0]
                tag = fastq.split("/")[-1].split(".")[0]
                lib0,lib1,bc,read = tag.split("_")
                key = lib0+"_"+lib1 + "_"+ bc
                if not key in self.fastqset.keys():
                    self.fastqset[key] = {}
                self.fastqset[key][read] = fastq
                if key not in self.fastqtag:
                    self.fastqtag.append(key)
            # check if all fastq are paired
            for k in self.fastqset.keys():
                if len(self.fastqset[k].keys()) != 2:
                    LG.error("fastq not paird for %s"%k)
            LG.info("loaded index for %d paird fastq"%len(self.fastqtag))
        pass

    def load_probe_annotation(self):
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
                tag = row[0] # we expect this to be unique name
                if tag not in seen:
                    seen.add(tag)
                else:
                    LG.error("Duplicated identifier:%s"%tag)
                    sys.exit()
                gene = row[1]
                region = row[11].split()
                chrm = region[0]
                start = int(region[1])
                end = int(region[2])
                fwseq = row[29]
                revseq = row[33]
                amplength = int(row[17])
                pb = Probe(chrm, start, end, tag, fwseq, revseq, amplength)
                pb.setGene(gene)
                self.probeset[pb.tag] = pb
        LG.info("Loaded %d probes" %(len(self.probeset.keys())))
    
    def build_hash(self):
        n = len(self.probeset.keys())
        assert n > 10, "probe information not loaded or too few probes"
        # this hash needs to be seq -> id
        seqdict = {}
        for k in self.probeset.keys():
            seqdict[self.probeset[k].fwseq]=k
            seqdict[self.probeset[k].revseq]=k
        self.primerHash = SeqHash(params["kmersize"], seqdict)
        LG.info("build hash table for %d probes using kmer size %d"%(n,params["kmersize"]))
        
    def test_hash(self):
        read = "NCACACATGCCCTCATTCTAGG" # note we introduced an erorr at the centre
        ids = self.primerHash.seedkmer(read)
        LG.info("expecting ID = 1")
        LG.info("found %s"%(str(ids)))
        
    def scan_fastq(self, index):
        tag = self.fastqtag[index]
        #tag = "WTCHG_98544_03"
        fastq_read1 = self.fastqset[tag]['1']
        fastq_read2 = self.fastqset[tag]['2']
        LG.info("read1 file %s "%fastq_read1)
        LG.info("read2 file %s "%fastq_read2)
        
        fastq_out_read1 = params["outdir"] + "/" + tag  + "_1.qc.gz"
        fastq_out_read2 = params["outdir"] + "/" + tag  + "_2.qc.gz"

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
        
        minn = params["minreadlength"]
        probe_search_length = params["probe_search_length"]
        start = params["trim_read"]

        counter = Counter(tag)
        ampcounter = {}
        mismatch_pairs = {}
        c = 0
        while True:
            c+=1
            if (c%1000 == 0):
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
            counter.nTotal_read += 1
            q1 = record1.letter_annotations["phred_quality"]
            q2 = record2.letter_annotations["phred_quality"]
            lowqual = False
            if ( sum(np.array(q1)< params["minQscore"]) > params["lowqualbase"]
                or sum(np.array(q2)< params["minQscore"]) > params["lowqualbase"]):
                counter.nLowQreads += 1
                lowqual = True
                continue
            n1 = len(record1.seq)
            n2 = len(record2.seq)
            if  n1!=n2 or n1 < minn or n2 < minn:
                counter.nShort_read += 1
                LG.debug("read too short or not consistent read1:%d read2:%d"%(n1,n2))
                continue
            id1 = self.primerHash.seedkmerMaxOnly(str(record1.seq[0:35]))
            id2 = self.primerHash.seedkmerMaxOnly(str(record2.seq[0:35]))
            if id1 == None or id2 == None:
                counter.nPrimer_NoMatch += 1
                continue
            if id1 == id2:
                counter.nQC_read += 1
                gene = self.probeset[id1].gene
                amplen = self.probeset[id1].amplength
                ## here is  to trim read two, which maynot be good as it makes it
                ## more difficult to align.
                # end2 = n2
                # LG.debug("id %s, n1 %d, n2 %d, amplen %d"%(id1,n1,n2,amplen))
                # if n1 + n2 > amplen:
                #     end2 = amplen - n1
                #     record2 = record2[start:end2]
                record1 = record1[start:]
                record2 = record2[start:]
                record1.id = id1 + ":" + gene + "|" + record1.id
                record2.id = id1 + ":" + gene + "|" + record2.id
                record1.description = ""
                record2.description = ""
                output_handle1.write(record1.format("fastq"))
                output_handle2.write(record2.format("fastq"))
                if not ampcounter.has_key(id1):
                    ampcounter[id1] = Counter(id1)
                ampcounter[id1].nTotal_read += 1
                if lowqual:
                    ampcounter[id1].nLowQreads += 1
            else:
                counter.nPrimer_MisMatch += 1
                if not mismatch_pairs.has_key(id1 + "-" + id2):
                    mismatch_pairs[id1 + "-" + id2] = 1
                else:
                    mismatch_pairs[id1 + "-" + id2] += 1

        R = {"global": counter, "amp": ampcounter, "mismatch": mismatch_pairs, "ampmeta":self.probeset}
        serial.cdmz(R, "%s/%s.picklez"%(params["outdir"], tag))

        handle1.close()
        handle2.close()
        output_handle1.close()
        output_handle2.close()

        LG.info("Completed scanning for %s"%tag)
        
    def summary_plots(self, index):
        tag = self.fastqtag[index]
        results = serial.clz("%s/%s.picklez"%(params["outdir"], tag))
        probes = results["ampmeta"]
        gl = results["global"]
        amp = results["amp"]
        mismatch_pairs = results["mismatch"]
        with open("%s/%s.csv"%(params["outdir"], tag), "wt") as outfh:
            gls = [gl.nTotal_read, gl.nQC_read,
                gl.nLowQreads, gl.nPrimer_MisMatch, gl.nPrimer_NoMatch]
            outfh.write(",".join([str(s) for s in gls]) + "\n")
        
            amps = []
            for a in  amp.keys():
                row = [probes[a].gene + "_" +a, amp[a].nTotal_read, amp[a].nLowQreads]
                amps.append(row)
                outfh.write(",".join([str(r) for r in row]) + "\n")
                
        from matplotlib.backends.backend_pdf import PdfPages
        figpath = "%s/%s.pdf"%(params["outdir"], tag)
        pdf_pages = PdfPages(figpath) # pdf pages
        
        fig = plt.figure()
        ax = plt.subplot(111)
        xpos = np.arange(1,len(gls)+1)
        xlabels = ["Total", "QC", "LowQual", "Primer_misMatch", "Primer_noMatch"]
        width = 0.5
        ax.bar(xpos, gls, width, color='k')
        ax.set_xticks(xpos+width/2)
        ax.set_xticklabels(xlabels, rotation='vertical', size='small')
        pdf_pages.savefig(fig,  bbox_inches='tight')
        
        # default scale is 1 in your original case, scales with other cases:
        widthscale = len(amps)/40
        figsize = (4*widthscale,6) # fig size in inches (width,height)
        fig = plt.figure(figsize = figsize)
#        fig = plt.figure()
        ax = plt.subplot(111)
        xpos = np.arange(1,len(amps)+1)
        labels = [i[0] for i in amps]
        atotals = [i[1] for i in amps]
        alowquals = [i[2] for i in amps]
        
        width = 0.3
        ax.bar(xpos, atotals, width, color='r', label="Total", edgecolor = "none")
        ax.bar(xpos+width, alowquals, width, color='k', label="LowQ", edgecolor = "none")
        ax.set_xticks(xpos)
        ax.set_xlim([0,len(xpos)+1])
        ax.set_xticklabels(labels, rotation='vertical', size=6)
        ax.legend()
        pdf_pages.savefig(fig,  bbox_inches='tight') # save the figure
    
        fig = plt.figure()
        ax = plt.subplot(211)
        ax.hist(atotals, 100, label="Total", color="r")
        ax.set_xlabel("Counts")
        ax.set_ylabel("Amplicons")
        ax.legend()
        ax = plt.subplot(212)
        ax.hist(alowquals, 100, label="Amplicons by low quality counts", color="k")
        ax.set_xlabel("Counts")
        ax.set_ylabel("Amplicons")
        ax.legend()
        pdf_pages.savefig(fig,  bbox_inches='tight') # save the figure
        
        pdf_pages.close() # close the file
        LG.info("Generated summary plots")
        plt.close()

    def globalsummary(self):
        LG.info("plotting global summary")
        nTotals = []
        nQCs = []
        nLowQreads = []
        nMisMatches = []
        nNoMatches = []
        for tag in self.fastqtag:
            results = serial.clz("%s/%s.picklez"%(params["outdir"], tag))
            gl = results["global"]
            nTotals.append(gl.nTotal_read)
            nQCs.append(gl.nQC_read)
            nLowQreads.append(gl.nLowQreads)
            nMisMatches.append(gl.nPrimer_MisMatch)
            nNoMatches.append(gl.nPrimer_NoMatch)

        fig = plt.figure()
        ax = plt.subplot(111)
        xpos = np.arange(len(self.fastqtag))
        xlabels = self.fastqtag
        width = 0.2
        scaler = 1
        ax.bar(xpos, nTotals, width, color='k', label="Total")
        ax.bar(xpos+scaler*width, nQCs, width, color='y', label="QC")
        ax.bar(xpos+scaler*2*width, nLowQreads, width, color='b', label="LowQ")
        ax.bar(xpos+scaler*3*width, nMisMatches, width, color='r', label="Mis")
        ax.bar(xpos+scaler*4*width, nNoMatches, width, color='g', label="NoMatch")
        ax.legend(prop={'size':7})
        ax.set_xticks(xpos+scaler*2*width)
        ax.set_xticklabels(xlabels, rotation='vertical', size='small')
        plt.savefig("%s/globalsummary.pdf"%(params["outdir"]),  bbox_inches='tight')
        LG.info("Generated global summary plot %s/globalsummary.pdf"%(params["outdir"]))
        
    def run_aligner(self, index):
        tag = self.fastqtag[index]
        lib0,lib1, bc = tag.split("_")
        lib = lib0 + "_" + lib1
        cmd = "/bin/bash " + params["aligner"] + " " + tag + " " + lib
        os.system(cmd)
        pass

class Pileup:
    
    # defined layout
    basepos = {"A":0,"C":1,"G":2,"T":3, "N":4}
    caselist = ['A', 'C', 'G', 'T', 'N', 'I', 'D']
    indel = {"I": 5, "D": 6}
    cases = {0:"A",1:"C",2:"G",3:"T", 4:"N", 5:"I", 6:"D"}
    nullidx = None
    n_reporting_cols = len(basepos.keys()) + len(indel.keys())
    readlength = params["readlength"]
    max_noncon_per_read = params["max_noncon_per_read"]
    
    def __init__(self):
        pass
    # pileup for each amplicon, only the target region
    # and the target reads
    # vcf like output
    # chrm, start, end, ref, alt, info, sample1 .... sampleN
    # apply caller
    
    def _baseMidx(self, seq):
        '''classify each base into ATGC'''
        idx = []
        for k in seq:
            b = k.upper()
            if self.basepos.has_key(b):
                idx.append(self.basepos[b])
            else:
                LG.error("found non ATGCN sequence in read %s"%seq)
                idx.append(self.nullidx)
        return idx

    def _cigarMidx(self, cigar):
        '''classify cigar operation, return index for M'''
        idx = []
        for c in cigar:
            if int(c) == 1:
                idx.append(self.indel["I"])
            elif int(c) == 2:
                idx.append(self.indel["D"])
            else:
                idx.append(self.nullidx)
        return idx
    
    def _cigarstr2seq(self, cigar):
        '''converte cigar string into an extended read of cigar letter'''
        cigar = re.findall(r'([0-9]+)([MIDNSHPX=])', cigar)
        cigar_read = ""
        for i in cigar:
            count, sym = i
            count = int(count)
            for j in range(count):
                cigar_read += sym
        return cigar_read

    def _cigartuple2seq(self, tuples):
        '''convert cigar tuple returned by pysam in tuple format, e.g
        (operation, counts) into a string representation.'''
        seq = ""
        for t in tuples:
            op, c = t
            seq += str(op)*c
        return seq
            
    def _fastqQuality(self, ascii_str):
        '''return the raw quality score Phred Q score'''
        return [ord(i)-33 for i in ascii_str]

    def countpileup(self):
        
        tag = "WTCHG_98544_01"
        bam = "/Users/zd1/mount/wimmhts/raindance/bams/%s/%s.sorted.bam"%(tag, tag)
#        outdir = "/Users/zd1/mount/wimmhts/raindance/bams/%s"%tag
        # pileup outdir
        outdir = "/Users/zd1/Work/projects2/raindance/pileup"
        samout = "%s/%s"%(params["outdir"], tag)
        
        #FGFR2_5
        amp = "FGFR2_5"
        ampid = "48"
        chrm= "chr10"
        
        ampstart = 123279459
        ampend = 123279658
        p0 = 123279492
        p1 = 123279643
                
        # amp = "FGFR2_5"
        # ampid = "49"
        # chrm= "chr10"
        # ampstart = 123279527
        # ampend = 123279702
        # p0 = 123279547
        # p1 = 123279684
        
        # amp = "FGFR2_5"
        # ampid = "50"
        # chrm= "chr10"
        # ampstart = 123279578
        # ampend = 123279734
        # p0 = 123279598
        # p1 = 123279714

        # amp specific dir, too many amps

        amptag = amp + '_' + ampid
        ampdir = "%s/amps/%s"%(samout, amptag)

        if not os.path.exists(ampdir):
            os.makedirs(ampdir)
        
                
        assert ampstart < p0 and ampend > p1
        targetsize = p1 - p0

        exp_read1_length = 100
        # keep the entire read1, trim read2
        exp_read2_length = 0
        if ampend - ampstart > 100:
            exp_read2_length = ampend - ampstart -100
        LG.debug("expected read2 length %d"%exp_read2_length)

        # these are relative positions to read1 and 2
        read1pos = np.array([str(i) for i in range(ampend - ampstart)])
        read2pos = np.array([str(i) for i in range(ampend - ampstart)])
        read1pos[100:] = "NA"
        read2pos[exp_read2_length:] = "NA"
        read2pos = read2pos[::-1]
        
        fa = Fasta(params["ref"])
        refseq =  fa.get_sequence_by_region(chrm, ampstart, ampend)
        refseq = refseq[:-1] # remove one bp, samtools gives inclusive for both ends
        LG.info("ref sequence %d basepairs"% len(refseq))
        LG.info("extracted reference sequence %s"%refseq)
        
        LG.info('running pileup for sample %s for amplicon %s(%s) at region %s:%d-%d' %(tag, amp, ampid, chrm, ampstart, ampend))
        outfile = "%s/%s.%s.MM.raw.out"%(ampdir, tag, ampid)
        
        #self._buld_consensus(tag, bam, outdir)
        tellyraw = self._count_pileup(tag, amp,  ampid, bam, chrm, ampstart, ampend, p0, p1, exp_read2_length, outfile)
        tellyraw.resolveMM()
        conseq = tellyraw.getConcensusSeq()
        conseqidx = tellyraw.getConcensusSeqIdx()
        LG.debug("concensus %d basepairs"%len(conseq))
        LG.debug(conseq)

        ndiff = 0
        for c in range(len(conseq)):
            if conseq[c] != refseq[c]:
                ndiff += 1

        if ndiff > 0.2*(ampend - ampstart):
            LG.error("Too many difference between the reference sequence and the consensus sequence")
            sys.exit(1)
            
        outfile = "%s/%s.%s.MM.qc.out"%(ampdir, tag, ampid)
        telly = self._count_pileup(tag, amp, ampid, bam, chrm, ampstart, ampend, p0, p1, outfile, exp_read2_length, conseqidx=conseqidx)
        telly.resolveMM()
        LG.debug(telly.nNonCon)
        LG.debug("hets:" + str(telly.idx_hets))
        
        # pass
        LG.info("total number of reads %d"%(telly.nTotal))
        LG.info("total number of reads after QC %d"%(telly.nQC))
        LG.info(" unpaired reads  %d"%(telly.nUnpaired))
        LG.info(" low mapping quality reads  %d"%(telly.nLowM))
        LG.info(" low quality bases normalised by reads %f"%(telly.nLowQbases*1.0/telly.nQC))
        LG.info(" reads filtered out due to having more than %d non consensus alleles  %d"%(params["max_noncon_per_read"],telly.nNonCon))
        
        LG.info(" found candidate sites  %d"%(len(telly.idx_candidates)))
        RM = np.hstack((np.array(telly.con_bases)[:,np.newaxis],
                        telly.coordinates[:, np.newaxis],
                        telly.MM,
                        np.array(telly.non_con_frq)[:,np.newaxis]))
        LG.info("finished pileup. exporting data ...")
        #serial.cdmz(telly, "%s/%s.telly.picklez"%(outdir, tag))
        telly.setRefseq(refseq)
        telly.setRead1pos(read1pos)
        telly.setRead2pos(read2pos)
        telly.store_data("%s/%s.MM.hd5"%(samout, tag))
        self._writeMM(RM, tag, amp, chrm, ampstart, ampend, p0, p1, refseq, exp_read1_length, exp_read2_length, read1pos, read2pos, "%s/%s.%s.MM.out"%(ampdir, tag, ampid))
        self._writeMM(RM, tag, amp, chrm, ampstart, ampend, p0, p1, refseq, exp_read1_length, exp_read2_length, read1pos, read2pos, "%s/%s.%s.MM.candidates.out"%(ampdir, tag, ampid), telly.idx_candidates)
        self._plotMM(tag, amp, ampid, telly, "%s/%s.%s.%s.pdf"%(ampdir, tag, amp, ampid))
        
        LG.info("Done.")
        
    def _align_chunk(self, readpos, readlen, ampstart, ampend):
        '''take the start and length of read, as well as
        the start and end of amplicion, work out the seq
        chunk of the read that overlaps the amplicon, and
        return indices for the chunk matrix, and indices
        for the read. If not overlapping, return None.'''
        readgpos = readpos + 1 # pos in pysam is 0-based
        
        #  ----t-----
        #------r--------
        if readgpos <= ampstart and readgpos + readlen >= ampend:
            rleft = 0 
            rright = ampend - ampstart
            seqstart = ampstart - readgpos
            seqend = ampend - readgpos
        #  ----t-----
        #------r---
        elif readgpos < ampstart and readgpos + readlen < ampend:
            rleft = 0
            rright = readgpos + readlen - ampstart
            seqstart = ampstart - readgpos
            seqend = readlen
        #  ----t-----
        #   ---r---
        elif readgpos >= ampstart and readgpos + readlen < ampend:
            rleft = readgpos - ampstart
            rright =  readgpos + readlen - ampstart
            seqstart = 0
            seqend = readlen
        #  ----t-----
        #   ---r------
        elif readgpos >= ampstart and readgpos + readlen >= ampend:
            rleft = readgpos - ampstart
            rright = ampend - ampstart
            seqstart = 0
            seqend = ampend - readgpos
        else:
            return None, None, None, None
        
        return rleft, rright, seqstart, seqend

    class Pileup_Telly:
        
        MM = None
        def __init__(self, chrm, ampstart, ampend, p0,p1, cases):
            self.ampstart = ampstart
            self.ampend = ampend
            self.nTotal= 0
            self.nQC = 0
            self.nUnpaired = 0
            self.nLowM = 0
            self.nLowQbases = 0
            self.nNonCon = 0
            self.coordinates = np.arange(ampstart, ampend)
            self.cases = cases
            self.chrm =  chrm
            self.p0 = p0
            self.p1 = p1
            self.con_bases = []
            self.con_base_colidx = []
            self.coverages = []
            self.non_con_frq = []
            self.idx_candidates = []
            self.idx_hets = []
            self.read12freq = []
            self.nQC_read1 = 0
            self.nQC_read2 = 0
        
        def setRefseq(self, refseq):
            self.refseq = refseq

        def setRead1pos(self, read1pos):
            self.read1pos = read1pos
            
        def setRead2pos(self, read2pos):
            self.read2pos = read2pos
        
        def getConcensusSeqIdx(self):
            return self.con_base_colidx
        
        def getConcensusSeq(self):
            return self.con_bases
        
        def resolveMM(self):
            for i in range(self.MM.shape[0]):
                row = self.MM[i,:]
                self.idx_con = np.argmax(row)
                self.con_base_colidx.append(self.idx_con)
                s = sum(row)
                nf = (s - row[self.idx_con])*1.0/s
                self.coverages.append(s)
                self.con_bases.append(self.cases[self.idx_con])
                self.non_con_frq.append(nf)
                if  nf > params["min_het_frq"]:
                    self.idx_hets.append(i)
                if  nf > params["min_non_consensus_freq"] and s > params["min_coverage"] and nf < params["min_het_frq"]:
                    self.idx_candidates.append(i)
            
        def store_data(self, samdatafile):
            '''here we store for each sample'''
            ampkey = "%s:%s-%s"%(self.chrm, self.p0, self.p1)
            try:
                fh = h5py.File(samdatafile, "w")
            
                dset = fh.create_dataset('%s/pileup'%ampkey, data=self.MM, chunks=True,compression='gzip')
                fh.create_dataset('%s/consensus_seq'%ampkey, data = self.con_bases, chunks=True,compression='gzip')
                fh.create_dataset('%s/noncon_frq'%ampkey, data = self.non_con_frq, chunks=True,compression='gzip')
                fh.create_dataset('%s/preliminary_can'%ampkey, data = self.idx_candidates, chunks=True,compression='gzip')
                fh.create_dataset('%s/hets'%ampkey, data = self.idx_hets, chunks=True,compression='gzip')
                fh.create_dataset('%s/read12frq'%ampkey, data = self.read12freq, chunks=True,compression='gzip')
                fh.create_dataset('%s/coverages'%ampkey, data = self.coverages, chunks=True,compression='gzip')
                fh.create_dataset('%s/refseq'%ampkey, data = list(self.refseq), chunks=True,compression='gzip')
                fh.create_dataset('%s/read1pos'%ampkey, data = self.read1pos, chunks=True,compression='gzip')
                fh.create_dataset('%s/read2pos'%ampkey, data = self.read2pos, chunks=True,compression='gzip')
            
                dset.attrs["nTotal"] = self.nTotal
                dset.attrs["nUnpaired"] = self.nUnpaired
                dset.attrs["nLowM"] = self.nLowM
                dset.attrs["nLowQbases"] = self.nLowQbases
                dset.attrs["nNonCon"] = self.nNonCon
                dset.attrs["nQC"] = self.nQC
                dset.attrs["nQC_read1"] = self.nQC_read1
                dset.attrs["nQC_read2"] = self.nQC_read2
            finally:
                fh.close()
            
    def _count_pileup(self, tag, amp,  ampid, bam, chrm, ampstart, ampend, p0, p1, exp_read2_length, outfile, conseqidx=None):
        
        samfile = pysam.Samfile("%s"%bam, "rb" )
        MM = np.zeros((ampend-ampstart, self.n_reporting_cols))
        MM_read1_only = np.zeros((ampend-ampstart, self.n_reporting_cols))
        MM_read2_only = np.zeros((ampend-ampstart, self.n_reporting_cols))
        
        M_read1 = np.zeros((self.readlength, self.n_reporting_cols))
        M_read2 = np.zeros((self.readlength, self.n_reporting_cols))
        read12freq = np.zeros((ampend-ampstart, 2))
        
        M_readbase=[]
        nTotal= 0
        nQC = 0
        nQC_read1 = 0
        nQC_read2 = 0
        nNonC_read1 = 0
        nNonC_read2 = 0
        nUnpaired = 0
        nLowM = 0
        nLowQbases = 0
        
        telly = self.Pileup_Telly(chrm, ampstart, ampend, p0, p1, self.cases)
        
        checkconsensus = False if conseqidx is None else True
        if checkconsensus:
            LG.info("checking concensus read")
            
        for alignedread in samfile.fetch(chrm, ampstart, ampend):
            
            if nTotal % 10000 == 0:
                LG.info("scanned %d reads"%nTotal)
            
            nTotal += 1
            rleft = None
            rright = None
            seq = None
            qual = None
            cigar = None
            
            # here we only look at the target amplicon
            readampid, readamptag = alignedread.qname.split("|")[0].split(":")
            if readampid != ampid or readamptag != amp:
                continue
            
            # only look at paired reads
            if not alignedread.is_paired:
                nUnpaired += 1
                continue
            # only look at reads with mapping quality above threshold
            if alignedread.mapq < params["mapQ_cutoff"]:
                nLowM += 1
                continue

            if alignedread.is_read2 and exp_read2_length == 0:
                continue
                
            # rleft, rright are respect to MM, seqstart/end are respect to read
            rleft, rright, seqstart, seqend = self._align_chunk(alignedread.pos, alignedread.rlen, ampstart, ampend)
            
            if rleft is None:
                continue
            
            seq = alignedread.seq[seqstart:seqend]
            qual =  self._fastqQuality(alignedread.qual[seqstart:seqend])
            cigar = self._cigartuple2seq(alignedread.cigar)[seqstart:seqend]
            
            idx = self._baseMidx(seq)
            idx_cigar = self._cigarMidx(cigar)
            
            # identify lowq base
            nbad = 0
            newseq = []
            for _b in range(len(qual)):
                if qual[_b] < params["baseQ_cutoff"]:
                    newseq.append("N")
                    nbad += 1
                else:
                    newseq.append(seq[_b])
            nLowQbases += nbad
            
            read12idx = 0 if alignedread.is_read1 else 1
            read12freq[rleft:rright, read12idx] += 1
            
            if checkconsensus:# remove reads with too many non con alleles
                n_non = 0
                for i in range(len(idx)):
                    if idx[i] != conseqidx[rleft+i]:
                        n_non += 1
                if n_non > self.max_noncon_per_read:
                    if alignedread.is_read1:
                        nNonC_read1 += 1
                    else:
                        nNonC_read2 += 1
                    telly.nNonCon += 1
                    continue

            if alignedread.is_read1:
                MM_read1_only = self._assignM(MM_read1_only, rleft, idx, idx_cigar)
            else:
                MM_read2_only = self._assignM(MM_read2_only, rleft, idx, idx_cigar)

            # here we consider how to combine read1/2
            offset = 0
            if alignedread.is_read2:
                # read2 needs to start after the end of the read1
                if alignedread.pos < alignedread.mpos + alignedread.rlen:
                    offset = alignedread.mpos + alignedread.rlen - alignedread.pos
                rleft += offset
                seqstart += offset

            idx = idx[offset:]
            idx_cigar = idx_cigar[offset:]
            
            for i in range(len(idx)):
                MM[rleft + i,idx[i]] += 1
                if idx_cigar[i] != self.nullidx:
                    MM[rleft+i,idx_cigar[i]] += 1
                    
            nQC += 1
            if alignedread.is_read1:
                nQC_read1 += 1
            else:
                nQC_read2 += 1
                
        samfile.close()
        LG.info("nQC_read1 %d nQC_read2 %d"%(nQC_read1, nQC_read2))
        LG.info("nNonC_read1 %d nNonC_read2 %d"%(nNonC_read1, nNonC_read2))
        
        telly.nTotal = nTotal
        telly.nQC = nQC
        telly.nLowQbases = nLowQbases
        telly.nLowM = nLowM
        telly.nUnpaired = nUnpaired
        telly.MM = MM
        telly.MM_read1_only = MM_read1_only
        telly.MM_read2_only = MM_read2_only
        telly.read12freq = read12freq
        telly.nQC_read1 = nQC_read1
        telly.nQC_read2 = nQC_read2
        
        return telly
    
    def _assignM(self, MM, rleft, idx, idx_cigar):
        for i in range(len(idx)):
            MM[rleft + i,idx[i]] += 1
            if idx_cigar[i] != self.nullidx:
                MM[rleft+i,idx_cigar[i]] += 1
        return MM
        
    def _plotMM(self, sam, amp, ampid, telly, imageout):
        chrm = telly.chrm
        ampstart = telly.ampstart
        ampend = telly.ampend
        p0 = telly.p0
        p1 = telly.p1
        con_base_colidx = telly.con_base_colidx
        con_bases = telly.con_bases
        coverages = telly.coverages
        read12freq = telly.read12freq
        idx_hets = telly.idx_hets
        coordinates = telly.coordinates
        
        LG.info("Generating summary plot")
        from matplotlib.backends.backend_pdf import PdfPages
        pdf_pages = PdfPages(imageout) # pdf pages

        allMM = [telly.MM, telly.MM_read1_only, telly.MM_read2_only]

        for MM in allMM:
            
            for i in range(MM.shape[0]):
                MM[i, con_base_colidx[i]] = 0 # don't show concensus base (number too large)
            MM[idx_hets,:] = 0 # don't show hets (number too large)
            fig = plt.figure() # plot on the first page
            ax = plt.subplot(111)
            xpos = np.arange(ampstart,ampend)
            
            colors = ['r', 'b', 'g', 'k', 'y', 'm', 'c']
            refcolormap = {self.caselist[i]:colors[i] for i in range(len(self.caselist))}
            
            for i in range(7):
                ax.plot(xpos, MM[:,i], marker=".", markersize=1, linestyle="None", color="w")
            for i in range(MM.shape[0]):
                for c in range(MM.shape[1]):
                    ax.annotate(self.cases[c], xy=(coordinates[i], MM[i,c]),
                            xytext=(coordinates[i], MM[i,c]),
                            xycoords='data', color=refcolormap[con_bases[i]], size=9)

            patches = []
            for c in self.caselist:
                patches.append(mpatches.Patch(color=refcolormap[c], label=c))
            ax.legend(handles=patches, prop={'size':7}, title="Concensus Allele")
            plt.setp(ax.get_legend().get_title(),fontsize='xx-small')
            
            ax.set_title("Sam %s, Amp %s ID %s (%s:%d-%d)"%(sam, amp, ampid, chrm, p0, p1),
                     fontdict={'fontsize': 10})
            ax.axvline(x=p0, color='k')
            ax.axvline(x=p1, color='k')
        
            camblue = (0.0627451 ,  0.39215686,  0.43921569)
            ax2 = ax.twinx()
            ax2.set_ylabel('Coverages', color= camblue)
            ax2.plot(xpos, coverages, linestyle= "-", color=camblue)
            ax2.set_ylim(bottom = 0, top = np.max(coverages)*1.1)
            ax2.set_xlim(right = ampend + 40)
            for tl in ax2.get_yticklabels():
                tl.set_color(camblue)
            
            pdf_pages.savefig(fig) # save the figure
            plt.clf()
            
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(xpos, read12freq[:,0], 'r-', label="read1")
        ax.plot(xpos, read12freq[:,1], 'y-', label="read2")
        ax.legend()
        # plot stuff
        pdf_pages.savefig(fig, bbox_inches='tight') # save the figure
        pdf_pages.close() # close the file
        LG.info("Done plotting")
        
    def _writeMM(self, MM, sam, amp, chrm, ampstart, ampend, p0, p1, refseq,
                 exp_read1_length, exp_read2_length, read1pos, read2pos, outfile, idx=None):
        
        ofh = file(outfile, "wb")
        header = ["Chrm", "Amp_Start","Amp_End", "Target_Start", 'Target_END',"Amplicon", "Sample", "Refseq",
                  "Concensus_Base", "Variant_Pos", "A", "C", "G", "T", "N", "I",
                  "D", "Non_Con_frq","Expected_Read1_Length", "Expected_Read2_length", "Pos_wrt_read1",  "Pos_wrt_read2"]
            
        ofh.write(",".join(header) + "\n")
        nrows = MM.shape[0]
        outrows = np.arange(nrows)
        if idx is not None:
            outrows = outrows[idx]
        for i in outrows:
            line = ",".join([chrm, str(ampstart), str(ampend), str(p0), str(p1), amp, sam]) + ","
            line += refseq[i] + ","
            line += ",".join([str(s) for s in MM[i,:]]) + ","
            line += str(exp_read1_length) + "," +  str(exp_read2_length) + ","
            line += str(read1pos[i]) + "," + str(read2pos[i])
            line += "\n"
            ofh.write(line)
        ofh.close()

        
class Fasta:
    # wrap sequence at chunck size
    chunksize = 80
    
    def __init__(self, ref):
        self.ref = ref
        assert os.path.exists(self.ref), "reference file not found"
        assert os.path.exists("%s.fai"%self.ref), "reference file not indexed"
        pass
    def get_sequence_by_region(self, chrm, start, end):
        '''use samtools faidx to retrieve ref data'''
        p1 = subprocess.Popen(["samtools","faidx",
                               self.ref,
                               "%s:%s-%s"%(chrm, start, end)],
                               stdout=subprocess.PIPE )
        seq = ""
        first = True
        while (p1.poll() == None):
            if p1.stdout == None:
                break
            for line in p1.stdout:
                if first:
                    first = False
                    continue
                line = line.strip()
                seq += line
            if not line:
                break
        return seq
    
class SummaryCounter:
    def __init__(self):
        self.nTotal = 0
        self.nProbeSequence = 0
        self.nProbeBothEnds = 0
        self.nExpectedLength = 0
        self.flen = [] # fragment length
    def summary2str(self, sep = ","):
        row = [self.nTotal, self.nProbeSequence, self.nProbeBothEnds,
                self.nExpectedLength]
        return sep.join([str(s) for s in row])
        
class SeqSample:
    def __init__(self, libname, barcode, read, path):
        self.libname = libname
        self.barcode = barcode
        self.read = read # this is 1 or 2
        self.path = path # file fastq file path
        pass
    
class Read:
    '''sequencing read of a PCR product'''
    def __init__(self, seq):
        self.seq = seq
    def getSeq(self):
        return self.seq
    def setSeq(self, seq):
        self.seq = seq
                
class Probe:
    '''we also call it amplicon'''
    def __init__(self, chrm, start, end,
                 tag, fwseq, revseq, amplength):
        self.chrm = chrm
        self.start = start
        self.end = end
        self.tag = tag
        self.fwseq = fwseq
        self.revseq = revseq
        self.amplength = amplength
        self.primerseq = None
        self.primerend = None
    def setGene(self, gene):
        self.gene = gene


class SeqSam:
    def __init__(self, flowcell, lane, samplebc, individual, slide, chunck):
        # chunk is "Testis Sample # in the sample meta data file"
        self.flowcell = flowcell
        self.lane = lane
        self.samplebc = samplebc
        self.individual = individual
        self.slide = slide
        self.chunck = chunck

                    
class SeqHash:
    '''build a kmer hash'''
    def __init__(self, kmersize, sequence_dict):
        self.kmersize = kmersize
        self._buildHash(sequence_dict)
        
    def _buildHash(self, sequence_dict):
        '''sequence_dict contains seq -> id pairs'''
        resulthash = {}
        for seq in sequence_dict.keys():
            seqlen = len(seq)
            if seqlen < self.kmersize:
                LG.error("can't use kmer that"
                         " is larger then the seuqence %s (kmersize%d)"%(seq, self.kmersize))
                sys.exit(1)
            for i in range(seqlen - self.kmersize):
                start = i
                end = i + self.kmersize
                kmer = seq[start:end]
                if resulthash.has_key(kmer):
                    resulthash[kmer].append(sequence_dict[seq])
                else:
                    resulthash[kmer] = [sequence_dict[seq]]
        self.hashtable = resulthash
        
    def _getKmers(self, read):
        '''break a read into kmers'''
        kmers = set()
        for i in range(len(read) - self.kmersize):
            start = i
            end = i + self.kmersize
            kmer = read[start:end]
            kmers.add(kmer)
        return kmers
        
    def seedkmer(self, read):
        kmers = self._getKmers(read)
        ids={}
        for kmer in kmers:
            if self.hashtable.has_key(kmer):
                for sam in self.hashtable[kmer]:
                    if not ids.has_key(sam):
                        ids[sam] = 1
                    else:
                        ids[sam] += 1
        return ids
    
    def seedkmerMaxOnly(self, read):
        '''consider the one with maximum amount of kmer match the match for
        the read'''
        allmatches = self.seedkmer(read)
        n_most = 0
        best=None
        if len(allmatches.keys()) > 0:
            for key in allmatches.keys():
                if allmatches[key] > n_most:
                    n_most = allmatches[key]
                    best = key
        return best
    
class Counter:
    nTotal_read = 0
    nShort_read = 0
    nQC_read = 0
    nPrimer_NoMatch = 0
    nPrimer_MisMatch = 0
    nLowQreads = 0
    def __init__(self, tag):
        self.tag = tag
    
# below are tests only
def testseed():
    "test read 45bp long"
    cigar = re.findall(r'([0-9]+)([MIDNSHPX=])', "1D66M34M")
    print cigar
    cigar_read = ""
    for i in cigar:
        count, sym = i
        count = int(count)
        for j in range(count):
            cigar_read += sym
        pass
    print cigar_read

def test():
    fa = Fasta("/Users/zd1/mount/wimmhts/common/ref/old/hg19_GRCh37/hg19.fa")
    print fa.get_sequence_by_region("chr10", "123279527", "123279702")
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Raindance analysis')
    parser.add_argument('--fastq', action='store_true', default=False)
    parser.add_argument('--pileup', action='store_true', default=False)
    parser.add_argument('--sample', required=True)
    args = parser.parse_args()
    print args
    
    # pr = cProfile.Profile()
    # pr.enable()
    
    if args.fastq:
        # runner = Analysis()
        # runner.load_seq_index()
        # runner.load_probe_annotation()
        # runner.build_hash()
        # runner.test_hash()
        # index = 0
        # runner.scan_fastq(index)
        # runner.summary_plots(index)
        print "run fastq"
    elif args.pileup:
        print "run pileup"
        # runner = Pileup()
        # runner.countpileup()
    else:
        print "need to specify --fastq or --pileup"
        sys.exit(1)
    # pr.disable()
    # s = StringIO.StringIO()
    # sortby = 'cumulative'
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print s.getvalue()
        

#    run()
