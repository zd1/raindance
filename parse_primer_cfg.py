
'''This is confirguration file
'''

# WTCHG
# param = {
#     "ref": "/well/mcvean/zd1/common/ref/old/hg19_GRCh37/hg19.fa",
#     "samplelist" : "/well/mcvean/zd1/raindance/meta/samplelist", # WT
#     #"samplelist" : "/Users/zd1/volumn/wt/raindance/meta/samplelist", # local
#     "fastqindex":"/well/mcvean/zd1/raindance/meta/fastq.index",
#     "table" : "/well/mcvean/zd1/raindance/meta/amplicons.csv", 
#     "ampshortlist": "/users/mcvean/zd1/workspace/raindance/amp_shortlist",
#     "fastq_outdir" : "/well/mcvean/zd1/raindance/fastq.qc",

#     # calling R script file
#     "call": "/users/mcvean/zd1/volumn/code/raindance/call.R",
#     #"call_out": "/users/mcvean/zd1/volumn/raindance/sum/calls_header.csv",
#     "call_out": "/Users/zd1/volumn/wt/raindance/sum/calls_header.csv",
#     #"ampkeys":"/Users/zd1/volumn/wt/raindance/pileups/all/keys", # local
#     "ampkeys":"/well/mcvean/zd1/raindance/meta/keys", # WT

#     # black lists
#     "samblack": "/users/mcvean/zd1/volumn/raindance/meta/sammut.blacklist",
#     "ampblack": "/users/mcvean/zd1/volumn/raindance/meta/amp.blacklist",
#     # "samblack": "/Users/zd1/cloud/data/raindance/meta/sammut.blacklist",
#     # "ampblack": "/Users/zd1/cloud/data/raindance/meta/amp.blacklist",

#     # filter
#     "p": 20,
#     "minncc": 10,
#     "minmut": 10,
#     #"vcfft": "/Users/zd1/cloud/data/raindance/pileup/sum/calls/all.flt.vcf.gz",
#     "vcfft": "/users/mcvean/zd1/volumn/raindance/sum/all.flt.vcf.gz",
#     "vcfanno":"/Users/zd1/volumn/wt/raindance/call/all.flt.vcf.anno.vcf.gz.vep",
    
#     # parameters:
#     "kmersize" : 10,
#     "lowqualbase" : 20, # reads with more than this amount of bad bases are considered low quality
#     "minQscore" : 20, # reads with more than this amount of bad bases are considered low quality
#     "readlength": 100,
    
#     # fastq filters
#     "minreadlength": 50,
#     "probe_search_length": 25,
#     "trim_read": 0,

#     # BAM directory. each sample is stored in bamdir/sample/sample.sorted.bam
#     "bamdir" : "/well/mcvean/zd1/raindance/bams",
    
#     # pileup out dirs
#     "pileup_outdir": "/well/mcvean/zd1/raindance/pileups",
#     "pileup_sum": "/users/mcvean/zd1/volumn/raindance/sum/all.hd5",
#     "pileup_agg": "/users/mcvean/zd1/volumn/raindance/sum/", # WT
#     "vcf": "/users/mcvean/zd1/volumn/raindance/sum//all.vcf.gz", #WT
#     # "pileup_agg": "/Users/zd1/cloud/data/raindance/pileup/sum", # local
#     # "vcf": "/Users/zd1/cloud/data/raindance/pileup/sum/all.vcf.gz", # local
    
#     # pileup filters
#     "baseQ_cutoff": 0, # use same cutoff as the one used to process fastq file
#     "mapQ_cutoff": 0, # 0-255, -10log10{mapping position is wrong}
#     "max_noncon_per_read": 10, # This is for identifying low quality reads. 
    
#     # pileup exporting
#     "min_coverage" : 5000, # min coverage
#     "min_non_consensus_freq": 0.005, # counts of all other bases together need to reach this frequency to be considered as elevated. 
#     "min_het_frq" : 0.20, # considered a het if non consensus counts exceed this threshold

#     #abet
#     "alphabet": ["A", "C", "G", "T"], 

#     #external data
#     "cosmic":"/users/mcvean/zd1/volumn/common/cosmic/cosmic.vcf.gz"
    
# }


# local
# WTCHG

param = {
    "ref": "/Users/zd1/volumn/wt/common/ref/old/hg19_GRCh37/hg19.fa",
    "samplelist" : "/Users/zd1/volumn/wt/raindance/meta/samplelist",
    "fastqindex":"/Users/zd1/volumn/wt/raindance/meta/fastq.index",
    "table" : "/Users/zd1/volumn/wt/raindance/meta/amplicons.csv", 
    "ampshortlist": "/users/mcvean/zd1/workspace/raindance/amp_shortlist",
    "fastq_outdir" : "/Users/zd1/volumn/wt/raindance/fastq.qc/",
    # parameters:
    "kmersize" : 10,
    "lowqualbase" : 20, # reads with more than this amount of bad bases are considered low quality
    "minQscore" : 20, # reads with more than this amount of bad bases are considered low quality
    "readlength": 100,
    
    # fastq filters
    "minreadlength": 50,
    "probe_search_length": 25,
    "trim_read": 0,

    # BAM directory. each sample is stored in bamdir/sample/sample.sorted.bam
    "bamdir" : "/Users/zd1/volumn/wimm/raindance/bams",
    
    # pileup out dirs
    # "pileup_outdir": "/Users/zd1/volumn/wimm/raindance/pileups",
    "pileup_outdir":"/Users/zd1/volumn/wt/raindance/pileups",
        
    # pileup filters
    "baseQ_cutoff": 0, # use same cutoff as the one used to process fastq file
    "mapQ_cutoff": 0, # 0-255, -10log10{mapping position is wrong}
    "max_noncon_per_read": 10, # This is for identifying low quality reads. 
    
    # pileup exporting
    "min_coverage" : 5000, # min coverage
    "min_non_consensus_freq": 0.005, # counts of all other bases together need to reach this frequency to be considered as elevated. 
    "min_het_frq" : 0.10, # considered a het if non consensus counts exceed this threshold
}
