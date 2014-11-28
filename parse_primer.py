import sys
import os
sys.path.append(os.environ['MYPYTHONPATH'])


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
    pa = Analysis()
    pa.load_probe_annotation()
    runner = Pileup()
    runner.pileupProbes(pa.probeset)
    
def run():
    run_pileup()
    
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
    fastqtag = [] # keys only in same order as input. e.g.  WTCHG_107218_01
    
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
                tag = row[0] # we expect this to be unique name, here we use the Amplicon ID
                if tag not in seen:
                    seen.add(tag)
                else:
                    LG.error("Duplicated identifier:%s"%tag)
                    sys.exit()
                gene = row[1]
                region = row[11].split() # this is the target region
                chrm = region[0]
                start = int(region[1])
                end = int(region[2])
                ampstart = int(row[22].replace(",",""))
                ampend = int(row[23].replace(",",""))
                fwseq = row[29]
                revseq = row[33]
                amplength = int(row[17])
                pb = Probe(chrm, start, end, ampstart, ampend, tag, fwseq, revseq, amplength)
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
        
    def scan_fastq(self, sample):
        #tag = self.fastqtag[index]
        #tag = "WTCHG_98544_03"
        tag = sample
        fastq_read1 = self.fastqset[tag]['1']
        fastq_read2 = self.fastqset[tag]['2']
        LG.info("read1 file %s "%fastq_read1)
        LG.info("read2 file %s "%fastq_read2)
        
        fastq_out_read1 = params["fastq_outdir"] + "/" + tag  + "_1.qc.gz"
        fastq_out_read2 = params["fastq_outdir"] + "/" + tag  + "_2.qc.gz"

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
        serial.cdmz(R, "%s/%s.picklez"%(params["fastq_outdir"], tag))

        handle1.close()
        handle2.close()
        output_handle1.close()
        output_handle2.close()

        LG.info("Completed scanning for %s"%tag)
        
    def summary_plots(self, sample):
        #tag = self.fastqtag[index]
        tag = sample
        results = serial.clz("%s/%s.picklez"%(params["fastq_outdir"], tag))
        probes = results["ampmeta"]
        gl = results["global"]
        amp = results["amp"]
        mismatch_pairs = results["mismatch"]
        with open("%s/%s.csv"%(params["fastq_outdir"], tag), "wt") as outfh:
            gls = [gl.nTotal_read, gl.nQC_read,
                gl.nLowQreads, gl.nPrimer_MisMatch, gl.nPrimer_NoMatch]
            outfh.write(",".join([str(s) for s in gls]) + "\n")
        
            amps = []
            for a in  amp.keys():
                row = [probes[a].gene + "_" +a, amp[a].nTotal_read, amp[a].nLowQreads]
                amps.append(row)
                outfh.write(",".join([str(r) for r in row]) + "\n")
                
        from matplotlib.backends.backend_pdf import PdfPages
        figpath = "%s/%s.pdf"%(params["fastq_outdir"], tag)
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
            results = serial.clz("%s/%s.picklez"%(params["fastq_outdir"], tag))
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
        plt.savefig("%s/globalsummary.pdf"%(params["fastq_outdir"]),  bbox_inches='tight')
        LG.info("Generated global summary plot %s/globalsummary.pdf"%(params["fastq_outdir"]))
        
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
    nullidx = -1
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
            if c == self.nullidx:
                idx.append(self.nullidx)
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

    def loadsamples(self):
        self.samples = []
        ofh = open(params["samplelist"], "r")
        for line in ofh:
            s = line.strip()
            self.samples.append(s)
        ofh.close()

    def _sub_sample(self, samplelist):
        subsamples = []
        self.loadsamples()
        ofh = open(samplelist, "r")
        for line in ofh:
            s = line.strip()
            if s in self.samples:
                subsamples.append(s)
            else:
                LG.error("sample %s doesn't exist in sample list %s"%(s, params["samplelist"]))
        ofh.close()
        return subsamples
    
    def _sub__amp(self, probes, amplist):
        subamps = []
        ofh = open(amplist, "r")
        for line in ofh:
            s = line.strip()
            if s in probes:
                subamps.append(s)
            else:
                LG.error("amplicon %s doesn't exist in amplicon list %s"%(s, params["table"]))
        ofh.close()
        return subamps
        
    def combine_pileup(self, probeset, amplist=None, samplelist=None):
        '''merge amps across samples into an integrated dataset'''
        amps = []
        samples = []
        if amplist is not None:
            amps = self._sub__amp(probeset.keys(), amplist)
        else:
            amps = probeset.keys()
        if samplelist is not None:
            samples = self._sub_sample(samplelist)
        else:
            samples = self.samples

        out = "%s/all"%(params["pileup_outdir"])
        if not  os.path.exists(out):
            os.mkdir(out)
            
        h5pyout = "%s/all.hd5"%(out)
        fh = h5py.File(h5pyout, "w")
        
        variables = ["nUnpaired", "nLowM", "nLowQbases", "nNonCon", "nQC", "nQC_read1", "nQC_read2"]
        
        M_variables = np.zeros((len(samples), len(amps), len(variables)))

        for i in range(len(samples)):
            sample = samples[i]
            samout = "%s/%s"%(params["pileup_outdir"], sample)
            samh5pyout = "%s/%s.MM.hd5"%(samout, sample)
            try:
                samfh = h5py.File(samh5pyout, "r")
            
                for a in range(len(amps)):
                    amp = amps[a]
                    ampkey = "%s:%s:%s:%s-%s"%(probeset[amp].tag, probeset[amp].gene, probeset[amp].chrm, probeset[amp].p0, probeset[amp].p1)
                
                    ampstart = probeset[amp].ampstart
                    ampend = probeset[amp].ampend
                    ampend += 1
                    
                    refseq = np.zeros(ampend-ampstart)
                    
                    if ampkey not in fh:
                        grp = fh.create_group(ampkey)
                        MM_amp = np.zeros((len(samples), ampend-ampstart, len(self.caselist)))
                        MM_amp_read1 = np.zeros((len(samples), ampend-ampstart, len(self.caselist)))
                        MM_amp_read2 = np.zeros((len(samples), ampend-ampstart, len(self.caselist)))
                                            
                        conseq = np.zeros((len(samples), ampend-ampstart))
                        read1pos = np.zeros((len(samples), ampend-ampstart))
                        read2pos = np.zeros((len(samples), ampend-ampstart))
                        qual1 = np.zeros((len(samples), ampend-ampstart)) # check if needs to be normalised
                        qual2 = np.zeros((len(samples), ampend-ampstart))
                        coverages = np.zeros((len(samples), ampend-ampstart))
                        # create init matrix for these
                        grp.create_dataset('%s/pileup'%ampkey, data = MM_amp, chunks=True,compression='gzip')
                        grp.create_dataset('%s/pileup_read1'%ampkey, data = MM_amp_read1, chunks=True,compression='gzip')
                        grp.create_dataset('%s/pileup_read2'%ampkey, data = MM_amp_read2, chunks=True,compression='gzip')
                        grp.create_dataset('%s/qual_read1'%ampkey, data = qual1, chunks=True,compression='gzip')
                        grp.create_dataset('%s/qual_read2'%ampkey, data = qual2, chunks=True,compression='gzip')
                        grp.create_dataset('%s/coverages'%ampkey, data = coverages, chunks=True,compression='gzip')
                        grp.create_dataset('%s/consensus_seq'%ampkey, data = conseq, chunks=True,compression='gzip')
                        # these only need to be set once
                        grp.create_dataset('%s/refseq'%ampkey, data = samfh['%s/refseq'%ampkey], chunks=True,compression='gzip')
                        grp.create_dataset('%s/target_index'%ampkey, data = samfh['%s/target_index'%ampkey], chunks=True,compression='gzip')

                    # sam specific data
                    fh["%s/pileup_read1"%ampkey][i,:,:] = samfh['%s/pileup_read1'%ampkey][:,:]
                    fh["%s/pileup_read2"%ampkey][i,:,:] = samfh['%s/pileup_read2'%ampkey][:,:]
                    fh["%s/qual_read1"%ampkey][i,:,:] = samfh['%s/qual_read1'%ampkey][:,:]
                    fh["%s/qual_read2"%ampkey][i,:,:] = samfh['%s/qual_read2'%ampkey][:,:]
                    fh["%s/coverages"%ampkey][i,:,:] = samfh['%s/coverages'%ampkey][:,:]
                    fh["%s/consensus_seq"%ampkey][i,:,:] = samfh['%s/coverages'%ampkey][:,:]
                    
                    # summary variables
                    for v in range(len(variables)):
                        M_variables[i,a,v] = fh['%s/pileup'%ampkey].attrs[variables[v]]
                            
            finally:
                samfh.close()

        fh.close()

    def pileupProbes(self, probeset, sample, amplistfile = None):
        bam = params["bamdir"]+"/%s/%s.sorted.bam"%(sample, sample)
        if not os.path.exists(bam):
            LG.error("can't find bam file %s"%bam)
            sys.exit(1)
        
        samout = "%s/%s"%(params["pileup_outdir"], sample)
        h5pyout = "%s/%s.MM.hd5"%(samout, sample)
        if os.path.exists(h5pyout):
            os.remove(h5pyout)

        amps = []
        if amplistfile is not None:
            ofh = open(amplistfile)
            amps = []
            for line in ofh:
                a = line.strip()
                if probeset.has_key(a):
                    amps.append(a)
                else:
                    LG.error("amp %s not found in amp list %s"%(a, params["table"]))

        if len(amps) > 0:
            for pb in amps:
                self.countpileup(probeset[pb], sample, bam, samout, h5pyout)
        else:
            for pb in probeset.keys():
                self.countpileup(probeset[pb], sample, bam, samout, h5pyout)
        
    def pileupProbes_test(self, probeset):
        samples = ["WTCHG_98544_05", "WTCHG_98544_03", "WTCHG_98544_05", "WTCHG_98544_08", "WTCHG_98544_09", "WTCHG_98544_13"]
        for sample in samples:
            bam = "/Users/zd1/volumn/wimm/raindance/bams/%s/%s.sorted.bam"%(sample, sample)
            samout = "%s/%s"%(params["pileup_outdir"], sample)
            h5pyout = "%s/%s.MM.hd5"%(samout, sample)
            if os.path.exists(h5pyout):
                os.remove(h5pyout)
            for pb in probeset.keys():
                if pb not in ["175"]:
                    continue
                self.countpileup(probeset[pb], sample, bam, samout, h5pyout)
            break
        
    def countpileup(self, probe, sample, bam, samout, h5pyout):
        
        tag = sample
        amp = probe.gene
        ampid = probe.tag
        chrm= probe.chrm
        ampstart = probe.ampstart
        ampend = probe.ampend
        p0 = probe.start
        p1 = probe.end

        # amp specific dir, too many amps
        amptag = amp + '_' + ampid
        ampdir = "%s/amps/%s"%(samout, amptag)
        ampend += 1 # including both end basepairs.  
        
        if not os.path.exists(ampdir):
            os.makedirs(ampdir)
        
        assert ampstart < p0 and ampend > p1
        targetsize = p1 - p0

        # mid point. when overlapping read1 stops here and read2
        # continues on. Note that this position is 1-based. 
        midpoint = int((ampstart + ampend)*1.0/2)
        
        exp_read1_length = midpoint - ampstart
        exp_read2_length = ampend - midpoint
        LG.debug("expected read1 length %d"%exp_read1_length)
        LG.debug("expected read2 length %d"%exp_read2_length)

        # these are relative positions to read1 and 2
        read1pos = np.array([str(i) for i in range(ampend - ampstart)])
        read2pos = np.array([str(i) for i in range(ampend - ampstart)])
        read1pos[exp_read1_length:] = "NA"
        read2pos[exp_read2_length:] = "NA"
        read2pos = read2pos[::-1]
        
        fa = Fasta(params["ref"])
        refseq =  fa.get_sequence_by_region(chrm, ampstart, ampend)
        refseq = refseq[:-1] # remove one bp, samtools gives inclusive for both ends
        LG.info("ref sequence %d basepairs"% len(refseq))
        LG.info("extracted reference sequence %s"%refseq)
        
        LG.info('running pileup for sample %s for amplicon %s(%s) at region %s:%d-%d' %(tag, amp, ampid, chrm, ampstart, ampend))
        outfile = "%s/%s.%s.MM.raw.out"%(ampdir, tag, ampid)

        # initial scan to build a blacklist, consisting of fragment names of fragments having 
        # elevated non ref alleles in either read of the read pair. 
        tellyraw = self._count_pileup(tag, amp,  ampid, bam, chrm, ampstart,
                                      ampend, p0, p1, midpoint, outfile, conseq=refseq)
        LG.debug("Identified low quality fragments %d"%(len(tellyraw.badfragments)))

        # result scan. this scan removes fragments from black list. 
        outfile = "%s/%s.%s.MM.qc.out"%(ampdir, tag, ampid)
        telly = self._count_pileup(tag, amp, ampid, bam, chrm, ampstart, ampend,
                                   p0, p1, midpoint, outfile, removebyname = tellyraw.badfragments)
        telly.resolveMM()
        LG.debug(telly.nNonCon)
        LG.debug("hets:" + str(telly.idx_hets))

        telly.ampid = ampid
        telly.ampname = amp
        
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
        telly.setRefseq(refseq)
        telly.setRead1pos(read1pos)
        telly.setRead2pos(read2pos)
        telly.store_data(h5pyout)
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
            self.Q_read1 = None
            self.Q_read2 = None
            self.badfragments = set()
            self.ampid = ""
            self.ampname = ""
            self.MM_read1_only = None
            self.MM_read2_only = None

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

        def buildMMFromReadpair(self):
            # here we use a quality rule. from data we think this is
            # the most sensible rule to use. 
            self.MM = np.zeros(self.MM_read1_only.shape)
            permswitch = False
            for i in range(self.MM_read1_only.shape[0]):
                if permswitch:
                    self.MM[i,:] = self.MM_read2_only[i,:]
                    continue
                if self.Q_read1[i] > self.Q_read2[i]:
                    self.MM[i,:] = self.MM_read1_only[i,:]
                else:
                    permswitch = True
                    self.MM[i,:] = self.MM_read2_only[i,:]
        
        def resolveMM(self):
            self.buildMMFromReadpair()
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
            ampkey = "%s:%s:%s:%s-%s"%(self.ampid, self.ampname, self.chrm, self.p0, self.p1)
            targetindex = []
            region = range(self.ampstart - self.ampend)
            for i in range(len(region)):
                if region[i] >= self.p0 and region[i] <= self.p1:
                    targetindex.append(i)
                
            try:
                fh = h5py.File(samdatafile, "a")
                if ampkey in fh:
                    del fh[ampkey]
                    
                dset = fh.create_dataset('%s/pileup'%ampkey, data=self.MM, chunks=True,compression='gzip')
                fh.create_dataset('%s/pileup_read1'%ampkey, data = self.MM_read1_only, chunks=True,compression='gzip')
                fh.create_dataset('%s/pileup_read2'%ampkey, data = self.MM_read2_only, chunks=True,compression='gzip')
                fh.create_dataset('%s/consensus_seq'%ampkey, data = self.con_bases, chunks=True,compression='gzip')
                fh.create_dataset('%s/noncon_frq'%ampkey, data = self.non_con_frq, chunks=True,compression='gzip')
                fh.create_dataset('%s/preliminary_can'%ampkey, data = self.idx_candidates, chunks=True,compression='gzip')
                fh.create_dataset('%s/read12frq'%ampkey, data = self.read12freq, chunks=True,compression='gzip')
                fh.create_dataset('%s/coverages'%ampkey, data = self.coverages, chunks=True,compression='gzip')
                fh.create_dataset('%s/refseq'%ampkey, data = list(self.refseq), chunks=True,compression='gzip')
                fh.create_dataset('%s/read1pos'%ampkey, data = self.read1pos, chunks=True,compression='gzip')
                fh.create_dataset('%s/read2pos'%ampkey, data = self.read2pos, chunks=True,compression='gzip')
                fh.create_dataset('%s/qual_read1'%ampkey, data = self.Q_read1, chunks=True,compression='gzip')
                fh.create_dataset('%s/qual_read2'%ampkey, data = self.Q_read2, chunks=True,compression='gzip')
                fh.create_dataset('%s/target_index'%ampkey, data = targetindex, chunks=True,compression='gzip')
                
                dset.attrs["ampstart"] = self.ampstart
                dset.attrs["ampend"] = self.ampend
                dset.attrs["targetstart"] = self.p0
                dset.attrs["targetend"] = self.p1
                
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
            
    def _count_pileup(self, tag, amp,  ampid, bam, chrm, ampstart, ampend,
                       p0, p1, midpoint, outfile, conseq=None, removebyname = None):

        # if the blacklisted readname is already constructed, do not check consensus again. 
        if removebyname is not None and len(removebyname) > 0:
            conseq = None
            
        checkconsensus = False
        if conseq is not None:
            LG.info("checking concensus read")
            checkconsensus = True
            
        samfile = pysam.Samfile("%s"%bam, "rb" )
        MM_read1_only = np.zeros((ampend-ampstart, self.n_reporting_cols))
        MM_read2_only = np.zeros((ampend-ampstart, self.n_reporting_cols))

        # base counts from each read
        M_read1 = np.zeros((self.readlength, self.n_reporting_cols))
        M_read2 = np.zeros((self.readlength, self.n_reporting_cols))
        # qaulity score from each read
        Q_read1 = np.zeros(ampend-ampstart)
        Q_read2 = np.zeros(ampend-ampstart)
        
        read12freq = np.zeros((ampend-ampstart, 2))
        
        M_readbase=[]
        nTotal= 0
        nQC = 0
        nQC_read1 = 0
        nQC_read2 = 0
        nNonC_read1 = 0
        nNonC_read2 = 0
        nNoMatch = 0
        nUnpaired = 0
        nLowM = 0
        nLowQbases = 0
        nLowQreads_read1 = 0
        nLowQreads_read2 = 0
        
        badfragments = set() # store the read names
        
        telly = self.Pileup_Telly(chrm, ampstart, ampend, p0, p1, self.cases)
        
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
                nNoMatch += 1
                continue
            
            # only look at paired reads
            if not alignedread.is_paired:
                nUnpaired += 1
                continue
            # only look at reads with mapping quality above threshold
            if alignedread.mapq < params["mapQ_cutoff"]:
                nLowM += 1
                continue
            
            # rleft, rright are respect to MM, seqstart/end are respect to read
            rleft, rright, seqstart, seqend = self._align_chunk(alignedread.pos, alignedread.qlen, ampstart, ampend)
            
            if rleft is None:
                continue
            
            seq = alignedread.query[seqstart:seqend]
            n_unknown = len([c for c in seq if c == 'N'])
            if n_unknown > 5:
                if alignedread.is_read1:
                    nLowQreads_read1 += 1
                else :
                    nLowQreads_read2 += 1
                continue
            
            # this is the base quality
            qual =  self._fastqQuality(alignedread.qqual[seqstart:seqend])
            if len(alignedread.cigar) > 0:
                cigar = self._cigartuple2seq(alignedread.cigar) # cigar for the entire read
                cigar = cigar[alignedread.qstart:alignedread.qend] # cigar for the mapped region only
                cigar = cigar[seqstart:seqend]  # cigar for the overlapped region
            else:
                cigar = [self.nullidx]*(seqend-seqstart)
                
            # store the quality scores
            if alignedread.is_read1:
                Q_read1[rleft:rright] += qual
            else:
                Q_read2[rleft:rright] += qual
                
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
            
            # remove flagged reads
            if removebyname is not None:
                if alignedread.qname in removebyname:
                    if alignedread.is_read1:
                        nNonC_read1 += 1
                    else:
                        nNonC_read2 += 1
                    telly.nNonCon += 1
                    continue
                
            # remove reads with too many non con alleles
            if checkconsensus:
                n_non = 0
                for i in range(len(idx)):
                    if self.cases[idx[i]] != conseq[rleft+i]:
                        n_non += 1
                if n_non > self.max_noncon_per_read:
                    if alignedread.is_read1:
                        nNonC_read1 += 1
                    else:
                        nNonC_read2 += 1
                    telly.nNonCon += 1
                    badfragments.add(alignedread.qname)
                    continue
                
            if alignedread.is_read1:
                MM_read1_only = self._assignM(MM_read1_only, rleft, idx, idx_cigar)
            else:
                MM_read2_only = self._assignM(MM_read2_only, rleft, idx, idx_cigar)
        
            nQC += 1
            if alignedread.is_read1:
                nQC_read1 += 1
            else:
                nQC_read2 += 1
                
        samfile.close()
        LG.info("nNomatch (in region not this amplicon) read %d"%(nNoMatch))
        LG.info("unPaired read %d"%(nUnpaired))
        LG.info("low mapping quality read %d"%(nLowM))
        LG.info("nQC_read1 %d nQC_read2 %d"%(nQC_read1, nQC_read2))
        LG.info("nNonC_read1 %d nNonC_read2 %d"%(nNonC_read1, nNonC_read2))
        LG.info("low quality reads. read1 %d read2 %d"%(nLowQreads_read1, nLowQreads_read2))
        
        telly.nTotal = nTotal
        telly.nQC = nQC
        telly.nLowQbases = nLowQbases
        telly.nLowM = nLowM
        telly.nUnpaired = nUnpaired
        telly.MM_read1_only = MM_read1_only
        telly.MM_read2_only = MM_read2_only
        telly.read12freq = read12freq
        telly.nQC_read1 = nQC_read1
        telly.nQC_read2 = nQC_read2
        telly.Q_read1 = Q_read1
        telly.Q_read2 = Q_read2
        telly.badfragments= badfragments
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

        avgQ1 = telly.Q_read1*1.0/(read12freq[:,0]+1)
        avgQ1[100:] = 0 # qual for position beyond read length is 0
        avgQ2 = telly.Q_read2*1.0/(read12freq[:,1]+1)
        avgQ2[0:-100] = 0 # qual for position beyond read length is 0
                
        LG.info("Generating summary plot")
        from matplotlib.backends.backend_pdf import PdfPages
        pdf_pages = PdfPages(imageout) # pdf pages

        allMM = [telly.MM, telly.MM_read1_only, telly.MM_read2_only]
        ymax = 0
        for MM in allMM:
            for i in range(MM.shape[0]):
                MM[i, con_base_colidx[i]] = 0 # don't show concensus base (number too large)
            MM[idx_hets,:] = 0 # don't show hets (number too large)
            thismax = np.max(MM)
            if thismax > ymax:
                ymax = thismax
        ymax = ymax *1.1
        
        # plot non consensus allele frequencies for combined, read1 and read2
        for MM in allMM:
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
            ax.set_ylim(top=ymax)
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
            ax2.set_xlim(left=ampstart-10, right = ampend + 40)
            for tl in ax2.get_yticklabels():
                tl.set_color(camblue)
            
            pdf_pages.savefig(fig) # save the figure
            plt.clf()

        # plot the base quality
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(xpos, avgQ1, 'r-', label="read1")
        ax.plot(xpos, avgQ2, 'y-', label="read2")
        ax.axvline(x=p0, color='k')
        ax.axvline(x=p1, color='k')
        ax.set_xlim(left=ampstart-10, right = ampend + 40)
        ax.legend()
        ax.set_title("FASTQ base quality")
        # plot stuff
        pdf_pages.savefig(fig, bbox_inches='tight') # save the figure
        plt.clf()
        
        # plot read1 and read2 area
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(xpos, read12freq[:,0], 'r-', label="read1")
        ax.plot(xpos, read12freq[:,1], 'y-', label="read2")
        ax.legend()
        # plot stuff
        pdf_pages.savefig(fig, bbox_inches='tight') # save the figure
        plt.clf()
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
    def __init__(self, chrm, start, end, ampstart, ampend,
                 tag, fwseq, revseq, amplength):
        self.chrm = chrm
        self.start = start # this is the target region
        self.end = end # this is the target region
        self.ampstart = ampstart
        self.ampend = ampend
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
    parser.add_argument('--amplist', required=False, default= None)
    args = parser.parse_args()
    LG.info(args)
    
    sample = args.sample
    
    if args.fastq:
        LG.info("processing fastq")
        runner = Analysis()
        runner.load_seq_index()
        runner.load_probe_annotation()
        runner.build_hash()
        runner.scan_fastq(sample)
        runner.summary_plots(sample)
    elif args.pileup:
        LG.info("doing pileup")
        pa = Analysis()
        pa.load_probe_annotation()
        amplistfile = args.amplist
        runner = Pileup()
        runner.pileupProbes(pa.probeset, sample, amplistfile)
    else:
        print "need to specify --fastq or --pileup"
        sys.exit(1)

    # run()
