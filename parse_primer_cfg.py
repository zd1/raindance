'''This is confirguration file
'''

# WTCHG
param = {
    "ref": "/well/mcvean/zd1/common/ref/old/hg19_GRCh37/hg19.fa",
    "samplelist" : "/well/mcvean/zd1/raindance/meta/samplelist",
    "fastqindex":"/well/mcvean/zd1/raindance/meta/fastq.index",
    "table" : "/well/mcvean/zd1/raindance/meta/amplicons.csv", 
    "ampshortlist": "/users/mcvean/zd1/workspace/raindance/amp_shortlist",
    "fastq_outdir" : "/well/mcvean/zd1/raindance/fastq.qc",
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
    "bamdir" : "/well/mcvean/zd1/raindance/bams",
    
    # pileup out dirs
    "pileup_outdir": "/well/mcvean/zd1/raindance/pileups",
    "pileup_sum": "/users/mcvean/zd1/volumn/raindance/sum/all.hd5",
    
    # pileup filters
    "baseQ_cutoff": 0, # use same cutoff as the one used to process fastq file
    "mapQ_cutoff": 0, # 0-255, -10log10{mapping position is wrong}
    "max_noncon_per_read": 10, # This is for identifying low quality reads. 
    
    # pileup exporting
    "min_coverage" : 5000, # min coverage
    "min_non_consensus_freq": 0.005, # counts of all other bases together need to reach this frequency to be considered as elevated. 
    "min_het_frq" : 0.10, # considered a het if non consensus counts exceed this threshold
}

# WIMM
# param = {
#     "ref": "/home/crangen/zding/projects/common/ref/old/hg19_GRCh37/hg19.fa",
#     "samplelist" : "/home/crangen/zding/projects/raindance/meta/samplelist",
#     "fastqindex":"/home/crangen/zding/projects/raindance/meta/fastq.index",
#     "table" : "/home/crangen/zding/projects/raindance/meta/amplicons.csv",
#     "ampshortlist": "/home/crangen/zding/workspace/raindance/amp_shortlist",
#     "fastq_outdir" : "/home/crangen/zding/projects/raindance/fastq.qc",
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
#     "bamdir" : "/home/crangen/zding/projects/raindance/bams",
    
#     # pileup out dirs
#     "pileup_outdir": "/home/crangen/zding/projects/raindance/pileups",
    
#     # pileup filters
#     "baseQ_cutoff": 0, # use same cutoff as the one used to process fastq file
#     "mapQ_cutoff": 0, # 0-255, -10log10{mapping position is wrong}
#     "max_noncon_per_read": 10, # This is for identifying low quality reads. 
    
#     # pileup exporting
#     "min_coverage" : 5000, # min coverage
#     "min_non_consensus_freq": 0.005, # counts of all other bases together need to reach this frequency to be considered as elevated. 
#     "min_het_frq" : 0.10, # considered a het if non consensus counts exceed this threshold

# }

# local
# WTCHG

# param = {
#     "ref": "/Users/zd1/volumn/wimm/common/ref/old/hg19_GRCh37/hg19.fa",
#     "samplelist" : "/Users/zd1/volumn/wimm/raindance/meta/samplelist",
#     "fastqindex":"/Users/zd1/volumn/wimm/raindance/meta/fastq.index",
#     "table" : "/Users/zd1/volumn/wimm/raindance/meta/amplicons.csv", 
#     "ampshortlist": "/users/mcvean/zd1/workspace/raindance/amp_shortlist",
#     "fastq_outdir" : "/Users/zd1/volumn/wimm/raindance/fastq.qc",
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
#     "bamdir" : "/Users/zd1/volumn/wimm/raindance/bams",
    
#     # pileup out dirs
#     # "pileup_outdir": "/Users/zd1/volumn/wimm/raindance/pileups",
#     "pileup_outdir":"/Users/zd1/volumn/wt/raindance/pileups",
        
#     # pileup filters
#     "baseQ_cutoff": 0, # use same cutoff as the one used to process fastq file
#     "mapQ_cutoff": 0, # 0-255, -10log10{mapping position is wrong}
#     "max_noncon_per_read": 10, # This is for identifying low quality reads. 
    
#     # pileup exporting
#     "min_coverage" : 5000, # min coverage
#     "min_non_consensus_freq": 0.005, # counts of all other bases together need to reach this frequency to be considered as elevated. 
#     "min_het_frq" : 0.10, # considered a het if non consensus counts exceed this threshold
# }
